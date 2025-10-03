#!/usr/bin/env bash
set -euo pipefail

# ===== 0) Project skeleton =====
ROOT=~/StMark/WGSeq/sa_polony
mkdir -p "${ROOT}"/{data,raw,qc,trimmed,assemblies,blast,amr,virulence,summary,logs,scripts,tmp}
cd "${ROOT}"

# ===== 1) Get reads (first 20 pairs from the provided list) =====
if [ ! -f SA_Polony_100_download.sh ]; then
  wget -q https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh
fi

# Make a smaller downloader for 20 samples, with resume/retry
sed -n '1,40p' SA_Polony_100_download.sh \
| sed -E 's#^curl #curl -C - --retry 3 #; s#^wget #wget -c --tries=3 #' \
> SA_Polony_20.sh

bash SA_Polony_20.sh

# Build sample lists and move fastqs into ./data
ls -1 *_1.fastq.gz | sed 's/_1\.fastq\.gz$//' | sort -u > samples_all.txt
head -n 20 samples_all.txt > samples_20.txt

# build a pairs table and move fastqs to data/
{
  echo -e "base\tR1\tR2"
  while read -r b; do
    echo -e "${b}\t${b}_1.fastq.gz\t${b}_2.fastq.gz"
  done < samples_20.txt
} > pairs_20.tsv

awk 'NR>1{print $2"\n"$3}' pairs_20.tsv | xargs -I{} mv -f {} data/

# ===== 2) QC + trimming (fastp) for the 20; we’ll analyze 10 =====
export THREADS=${THREADS:-8}

while read -r S; do
  b=$(basename "${S}")
  if [ -s trimmed/${b}_1.trim.fastq.gz ] && [ -s trimmed/${b}_2.trim.fastq.gz ]; then
    echo "[trim] Skip ${b} (found trimmed)"
    continue
  fi
  fastp \
    -i data/${b}_1.fastq.gz -I data/${b}_2.fastq.gz \
    -o trimmed/${b}_1.trim.fastq.gz -O trimmed/${b}_2.trim.fastq.gz \
    -w ${THREADS} --detect_adapter_for_pe \
    --html qc/${b}.fastp.html --json qc/${b}.fastp.json
done < samples_20.txt

head -n 10 samples_20.txt > samples_10.txt

# ===== 3) Assembly (SPAdes) for the 10 =====
mkdir -p assemblies logs
while read -r b; do
  if [ -s assemblies/${b}.contigs.fasta ]; then
    echo "[spades] Skip ${b} (contigs exist)"; continue
  fi
  echo "[spades] Assembling ${b} ..."
  spades.py \
    -1 trimmed/${b}_1.trim.fastq.gz \
    -2 trimmed/${b}_2.trim.fastq.gz \
    -o assemblies/${b} \
    -t ${THREADS} -m 16 --only-assembler \
    2>&1 | tee logs/spades_${b}.log
  cp assemblies/${b}/contigs.fasta assemblies/${b}.contigs.fasta
done < samples_10.txt

# ===== 4) Assembly sanity stats (no seqkit/quast needed) =====
python3 - <<'PY'
import os, glob, gzip
def readfa(p):
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
        if name: yield name, ''.join(seq)

os.makedirs('summary', exist_ok=True)
with open('summary/assembly_basic.tsv','w') as w:
    w.write("SAMPLE\tN_CONTIGS\tTOTAL_BP\n")
    for fa in sorted(glob.glob('assemblies/*.contigs.fasta')):
        s = os.path.basename(fa).replace('.contigs.fasta','')
        n, total = 0, 0
        for _,seq in readfa(fa):
            n += 1
            total += len(seq)
        w.write(f"{s}\t{n}\t{total}\n")
print("Wrote summary/assembly_basic.tsv")
PY

# ===== 5) Build BLAST query: largest contig per sample =====
mkdir -p blast
python3 - <<'PY'
import glob, gzip
def readfa(p):
    op = gzip.open if p.endswith('.gz') else open
    with op(p,'rt') as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name, seq = line[1:].strip().split()[0], []
            else: seq.append(line.strip())
        if name: yield name, ''.join(seq)

with open('blast/largest_contigs.fasta','w') as w:
    for fa in sorted(glob.glob('assemblies/*.contigs.fasta')):
        sample = fa.split('/')[-1].replace('.contigs.fasta','')
        best = max(readfa(fa), key=lambda x: len(x[1]))
        w.write(f">{sample}|{best[0]}\n{best[1]}\n")
print("Wrote blast/largest_contigs.fasta")
PY

# ===== 6) Remote BLAST (lean, genus-filtered) with nohup =====
# This runs in background so the session can drop without killing it.
if [ ! -s blast/blast_refseq.tsv ]; then
  nohup bash -lc '
    blastn -task megablast \
      -query blast/largest_contigs.fasta \
      -db refseq_genomic -remote \
      -entrez_query "txid1637[Organism:exp]" \
      -qcov_hsp_perc 80 -perc_identity 90 \
      -max_hsps 1 -max_target_seqs 3 -evalue 1e-40 \
      -outfmt "6 qseqid sacc pident length evalue bitscore staxids stitle" \
      -out blast/blast_refseq.tsv
  ' > logs/blast_refseq.stdout 2> logs/blast_refseq.err &
  echo $! > logs/blast_refseq.pid
  echo "[blast] Started background BLAST PID $(cat logs/blast_refseq.pid)"
else
  echo "[blast] Found existing blast/blast_refseq.tsv"
fi

# Reduce to one top hit per sample (safe to rerun anytime)
if [ -s blast/blast_refseq.tsv ]; then
  awk 'BEGIN{FS=OFS="\t"} {split($1,a,"|"); s=a[1]; if(!(s in seen)){print s,$8,$3,$4,$6,$2,$5; seen[s]=1}}' \
    blast/blast_refseq.tsv | sort -k1,1 \
    > blast/blast_top_hit.tsv
fi

# ===== 7) AMR (ABRicate / ResFinder) =====
# Prefer preinstalled central DB dir; else download to $HOME and set ABRICATE_DB
if abricate --list 2>/dev/null | awk 'NR>1{print $1}' | grep -qx resfinder; then
  export ABRICATE_DB=/opt/conda/db
else
  export ABRICATE_DB=~/db/abricate
  mkdir -p "$ABRICATE_DB"
  # download minimal DBs we need
  for db in resfinder vfdb; do
    abricate-get_db --dbdir "$ABRICATE_DB" --db "$db" || true
  done
  abricate --setupdb
fi

# Run ResFinder
abricate --datadir "$ABRICATE_DB" --db resfinder assemblies/*.contigs.fasta \
  > amr/resfinder_hits.tsv
abricate --summary amr/resfinder_hits.tsv > amr/resfinder_summary.tsv

# Quick prevalence table (assumes 10 isolates here)
awk 'NR>1{g[$6]++} END{
  printf "GENE\tISOLATES\tPREVALENCE(%%)\n";
  for(k in g) printf "%s\t%d\t%.1f\n",k,g[k],(g[k]/10)*100
}' amr/resfinder_hits.tsv | sort -k2,2nr \
> amr/resfinder_gene_prevalence.tsv

# ===== 8) Virulence / VFDB (incl. hly, plcA, plcB) =====
abricate --datadir "$ABRICATE_DB" --db vfdb assemblies/*.contigs.fasta \
  > virulence/vfdb_hits.tsv
abricate --summary virulence/vfdb_hits.tsv > virulence/vfdb_summary.tsv

awk 'BEGIN{FS=OFS="\t"} NR==1{next} $6 ~ /(^|,)(hly|plcA|plcB)(,|$)/ {print $1,$6,$10,$11}' \
  virulence/vfdb_hits.tsv \
| sort -k1,1 -k2,2 \
| awk 'BEGIN{print "SAMPLE\tGENE\t%COV\t%IDENTITY"}{print}' \
> virulence/toxin_hly_plcA_plcB.tsv

# ===== 9) Helpful “how is BLAST doing?” tips (print once) =====
echo
echo "BLAST status:"
if [ -f logs/blast_refseq.pid ]; then
  ps -fp "$(cat logs/blast_refseq.pid 2>/dev/null)" || echo "No running BLAST process"
fi
[ -f blast/blast_refseq.tsv ] && wc -l blast/blast_refseq.tsv || true
echo "Check: tail -n 50 logs/blast_refseq.err"
echo
echo "DONE (pipeline finished). If BLAST still running, it will keep writing to blast/blast_refseq.tsv."
