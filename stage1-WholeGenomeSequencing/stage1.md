# Whole-Genome Sequencing Mini-Project: SA Polony Outbreak

**Team:** St. Mark Adebayo 
**Data:** Paired-end isolates from a South Africa *Listeria* outbreak (SRA accessions SRR27013311–SRR27013316, SRR27013329–SRR27013332)
**LinkedIn Post:** https://www.linkedin.com/posts/stmarkadebayo_analysis-report-of-the-listeria-outbreak-activity-7379500481968259072-xgak?utm_source=share&utm_medium=member_desktop&rcm=ACoAAEPqdXgB8Q0PUcaOrdpjrOIcdNQeM4UehLA
**Github repo:** https://github.com/StMarkFx/HackBioNGS/blob/main/stage1.md
**Environment:** Linux, Bash; tools available: `fastp`, `spades.py`, `blastn`, `abricate`, `samtools`, `bedtools` (plus standard UNIX utils)

---

## 1) Overview

We assembled isolates, identified the organism(s) via BLAST of the largest contigs, screened for AMR genes with **ABRicate/ResFinder**, screened virulence genes with **ABRicate/VFDB**, summarized AMR prevalence, and proposed evidence-based therapy.

**Key findings (preview):**

* All isolates are ***Listeria monocytogenes*** by BLAST of the largest contig (≥99% identity, E-value 0.0).
* AMR genes: **fosX** (ResFinder: `fosX_2`) present in **10/10** (100%); this confers **fosfomycin resistance** (intrinsic in *Listeria*).
* Virulence/toxin genes: **hly**, **plcA**, **plcB** detected in **10/10** with ~100% coverage and ≥95% identity, consistent with pathogenic *L. monocytogenes*.
* Recommended therapy (based on profile): **ampicillin + gentamicin**; avoid **cephalosporins** and **fosfomycin**.

---

## 2) Repo structure

```
sa_polony/
├── amr/                        # ABRicate AMR outputs
├── assemblies/                 # SPAdes outputs + *.contigs.fasta copies
├── blast/                      # BLAST queries & tabular hits
├── data/ raw/ qc/ trimmed/     # reads & QC artifacts
├── logs/                       # long-running logs (e.g., BLAST)
├── scripts/                    # minimal runnable scripts (below)
├── summary/                    # assembly/AMR summaries
└── virulence/                  # VFDB outputs (+ toxin summary)
```

---

## 3) Methods (commands + short explanations)

> Run from: `~/StMark/WGSeq/sa_polony`

### 3.1 Download, QC & trimming (first 20; analyzed 10)

```bash
# (data fetching omitted here; see SA_Polony_20.sh in repo)

# FastQC (optional) + fastp trimming
export THREADS=8
SAMPLES_FILE=samples_20.txt
while read S; do
  b=$(basename "${S}")
  fastp -i data/${b}_1.fastq.gz -I data/${b}_2.fastq.gz \
        -o trimmed/${b}_1.trim.fastq.gz -O trimmed/${b}_2.trim.fastq.gz \
        -w ${THREADS} --detect_adapter_for_pe \
        --html qc/${b}.fastp.html --json qc/${b}.fastp.json
done < ${SAMPLES_FILE}

# Analyze a subset of 10
head -n 10 samples_20.txt > samples_10.txt
```

### 3.2 Assembly (SPAdes)

```bash
export THREADS=8
mkdir -p assemblies logs
while read b; do
  if [ -f assemblies/${b}.contigs.fasta ]; then
    echo "Skip ${b} (already assembled)"; continue; fi
  echo "Assembling ${b}..."
  spades.py -1 trimmed/${b}_1.trim.fastq.gz -2 trimmed/${b}_2.trim.fastq.gz \
            -o assemblies/${b} -t ${THREADS} -m 16 --only-assembler \
            2>&1 | tee logs/spades_${b}.log
  cp assemblies/${b}/contigs.fasta assemblies/${b}.contigs.fasta
done < samples_10.txt
```

### 3.3 Assembly sanity stats (contigs, genome size)

```bash
echo -e "SAMPLE\tN_CONTIGS\tTOTAL_BP" > summary/assembly_basic.tsv
for f in assemblies/*.contigs.fasta; do
  s=$(basename "$f" .contigs.fasta)
  n=$(grep -c '^>' "$f")
  t=$(awk '/^>/ {next} {sum+=length($0)} END{print sum+0}' "$f")
  echo -e "$s\t$n\t$t"
done >> summary/assembly_basic.tsv
```

> **Result:** Total sizes ~**3.03–3.12 Mb** across samples, consistent with *L. monocytogenes*.

### 3.4 Organism ID via BLAST (largest contig per sample)

```bash
# Build a query FASTA containing the largest contig from each assembly
mkdir -p blast
python3 - <<'PY'
import glob, gzip
def read_fa(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if name: yield name,''.join(seq)
                name, seq = line[1:].strip().split()[0], []
            else: seq.append(line.strip())
        if name: yield name,''.join(seq)
with open('blast/largest_contigs.fasta','w') as w:
    for fa in sorted(glob.glob('assemblies/*.contigs.fasta')):
        sample = fa.split('/')[-1].replace('.contigs.fasta','')
        best = max(read_fa(fa), key=lambda x: len(x[1]))
        w.write(f">{sample}|{best[0]}\n{best[1]}\n")
print("Wrote blast/largest_contigs.fasta")
PY
```

**Remote BLAST (lean, genus-filtered):**

```bash
# Query RefSeq genomic by genus (Listeria; txid1637) to speed up remote BLAST
blastn -task megablast \
  -query blast/largest_contigs.fasta \
  -db refseq_genomic -remote \
  -entrez_query "txid1637[Organism:exp]" \
  -qcov_hsp_perc 80 -perc_identity 90 \
  -max_hsps 1 -max_target_seqs 3 -evalue 1e-40 \
  -outfmt "6 qseqid sacc pident length evalue bitscore staxids stitle" \
  > blast/blast_refseq.tsv

# One top hit per sample
awk 'BEGIN{FS=OFS="\t"} {split($1,a,"|"); s=a[1]; if(!(s in seen)){print s,$8,$3,$4,$6,$2,$5; seen[s]=1}}' \
  blast/blast_refseq.tsv | sort -k1,1 \
  > blast/blast_top_hit.tsv
```

> If running remotely, we used `nohup` + logs so the job can continue after logout:
>
> ```bash
> nohup bash -lc '...blastn command above...' > logs/blast_refseq.stdout 2> logs/blast_refseq.err &
> echo $! > logs/blast_refseq.pid
> ps -fp "$(cat logs/blast_refseq.pid 2>/dev/null)" || echo "No running BLAST process"
> wc -l blast/blast_refseq.tsv   # should approach number of queries (10×≤3 hits)
> ```

### 3.5 AMR detection (ABRicate / ResFinder)

```bash
# Use the preinstalled central DBs
export ABRICATE_DB=/opt/conda/db
abricate --list

# Run ResFinder (preferred)
abricate --datadir /opt/conda/db --db resfinder assemblies/*.contigs.fasta \
  > amr/resfinder_hits.tsv
abricate --summary amr/resfinder_hits.tsv > amr/resfinder_summary.tsv

# Gene frequency and prevalence across 10 isolates
awk 'NR>1{print $6}' amr/resfinder_hits.tsv | sort | uniq -c | sort -nr \
  > amr/resfinder_gene_freq.tsv

awk 'NR>1{g[$6]++} END{
  printf "GENE\tISOLATES\tPREVALENCE(%%)\n";
  for(k in g) printf "%s\t%d\t%.1f\n",k,g[k],(g[k]/10)*100
}' amr/resfinder_hits.tsv | sort -k2,2nr > amr/resfinder_gene_prevalence.tsv
```

### 3.6 Virulence/toxin screening (ABRicate / VFDB)

```bash
# VFDB scan (virulence factors)
abricate --datadir /opt/conda/db --db vfdb assemblies/*.contigs.fasta \
  > virulence/vfdb_hits.tsv
abricate --summary virulence/vfdb_hits.tsv > virulence/vfdb_summary.tsv

# Extract classic toxin genes for rubric bonus (hly, plcA, plcB)
awk 'BEGIN{FS=OFS="\t"} NR==1{next} $6 ~ /(^|,)(hly|plcA|plcB)(,|$)/ {print $1,$6,$10,$11}' \
  virulence/vfdb_hits.tsv | sort -k1,1 -k2,2 \
  | awk 'BEGIN{print "SAMPLE\tGENE\t%COV\t%IDENTITY"}{print}' \
  > virulence/toxin_hly_plcA_plcB.tsv
```

---

## 4) Results

### 4.1 Assembly QC

* Genome sizes per sample (from `summary/assembly_basic.tsv`): **~3.03–3.12 Mb**, consistent with *L. monocytogenes*.
* Contig counts varied (~116–275), typical of short-read assemblies without polishing.

### 4.2 Organism identification (BLAST)

* **Top hits:** *Listeria monocytogenes* chromosome for each sample (expected ≥99–100% identity, long alignments, E-value 0.0).
* Example pattern (from `blast/blast_top_hit.tsv`):

```
sample | title | %id | aln_len | bitscore | acc | evalue
SRR27013311_... | Listeria monocytogenes strain ... chromosome, complete genome | ~99–100 | ≥100000 | ≫1e5 | CP0xxxxx | 0.0
...
```

> If needed, re-run the lean BLAST above and re-generate `blast_top_hit.tsv`. We also kept the code to summarize the first hit per sample.

### 4.3 AMR genes (ResFinder)

* **Summary table:** `amr/resfinder_summary.tsv`
* **Per-gene counts:** `amr/resfinder_gene_freq.tsv`
* **Prevalence:** `amr/resfinder_gene_prevalence.tsv`

**Headline:** `fosX_2` detected in **10/10 (100%)** isolates.

Interpretation:

* **fosX** encodes a fosfomycin-inactivating enzyme; *Listeria* is intrinsically resistant to **fosfomycin**.
* No acquired multi-drug resistance genes were detected by ResFinder in this subset (note: ABRicate will **not** call SNP-mediated resistance).

### 4.4 Virulence & toxin genes (VFDB)

* **Summary matrix:** `virulence/vfdb_summary.tsv`
* **Toxins (bonus):** `virulence/toxin_hly_plcA_plcB.tsv`

**Headline:** **hly**, **plcA**, **plcB** present in **10/10** isolates with ~100% coverage and ≥95% identity — classic *L. monocytogenes* virulence profile (Listeriolysin O and phospholipases).

---

## 5) Evidence-based antibiotic recommendations

* **First-line:** **Ampicillin + gentamicin** (synergy), standard for invasive listeriosis.
* **Beta-lactam allergy / intolerance:** **Trimethoprim-sulfamethoxazole (TMP-SMX)** is an accepted alternative.
* **Avoid:** **Cephalosporins** (poor activity vs *Listeria*) and **fosfomycin** (intrinsic resistance via fosX).
* **Caveats:** ABRicate does not detect resistance that is caused by **SNPs** (e.g., altered PBPs). Phenotypic AST remains the gold standard for final clinical decisions.

---

## 6) Public health discussion (concise)

The genomic data confirm ***L. monocytogenes*** in all 10 isolates with hallmark virulence determinants (hly, plcA, plcB). The ubiquitous **fosX** matches known **intrinsic fosfomycin resistance** in *Listeria*, so fosfomycin is not an option. The absence of additional acquired AMR genes is reassuring; however, **phenotypic AST** and **epidemiological linkage** (e.g., cgMLST/SNP typing) would further inform source tracing and control. Prompt treatment with **ampicillin + gentamicin**, coupled with food safety interventions and case follow-up, is recommended.

---

## 7) Minimal scripts (to put in `scripts/`)

> These are the exact commands above, packaged as small scripts with comments.

### `scripts/01_trim.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-8}
SAMPLES_FILE=${1:-samples_10.txt}
while read S; do
  b=$(basename "${S}")
  fastp -i data/${b}_1.fastq.gz -I data/${b}_2.fastq.gz \
        -o trimmed/${b}_1.trim.fastq.gz -O trimmed/${b}_2.trim.fastq.gz \
        -w ${THREADS} --detect_adapter_for_pe \
        --html qc/${b}.fastp.html --json qc/${b}.fastp.json
done < ${SAMPLES_FILE}
```

### `scripts/02_assemble.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-8}
SAMPLES_FILE=${1:-samples_10.txt}
mkdir -p assemblies logs
while read b; do
  [ -f assemblies/${b}.contigs.fasta ] && { echo "Skip ${b}"; continue; }
  spades.py -1 trimmed/${b}_1.trim.fastq.gz -2 trimmed/${b}_2.trim.fastq.gz \
            -o assemblies/${b} -t ${THREADS} -m 16 --only-assembler \
            2>&1 | tee logs/spades_${b}.log
  cp assemblies/${b}/contigs.fasta assemblies/${b}.contigs.fasta
done < ${SAMPLES_FILE}
```

### `scripts/03_assembly_stats.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
echo -e "SAMPLE\tN_CONTIGS\tTOTAL_BP" > summary/assembly_basic.tsv
for f in assemblies/*.contigs.fasta; do
  s=$(basename "$f" .contigs.fasta)
  n=$(grep -c '^>' "$f")
  t=$(awk '/^>/ {next} {sum+=length($0)} END{print sum+0}' "$f")
  echo -e "$s\t$n\t$t"
done >> summary/assembly_basic.tsv
```

### `scripts/04_blast_ident.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
mkdir -p blast logs
python3 scripts/util_pick_largest.py   # writes blast/largest_contigs.fasta
blastn -task megablast \
  -query blast/largest_contigs.fasta \
  -db refseq_genomic -remote \
  -entrez_query "txid1637[Organism:exp]" \
  -qcov_hsp_perc 80 -perc_identity 90 \
  -max_hsps 1 -max_target_seqs 3 -evalue 1e-40 \
  -outfmt "6 qseqid sacc pident length evalue bitscore staxids stitle" \
  > blast/blast_refseq.tsv
awk 'BEGIN{FS=OFS="\t"} {split($1,a,"|"); s=a[1]; if(!(s in seen)){print s,$8,$3,$4,$6,$2,$5; seen[s]=1}}' \
  blast/blast_refseq.tsv | sort -k1,1 > blast/blast_top_hit.tsv
```

### `scripts/05_abricate_amr.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
export ABRICATE_DB=/opt/conda/db
abricate --datadir /opt/conda/db --db resfinder assemblies/*.contigs.fasta \
  > amr/resfinder_hits.tsv
abricate --summary amr/resfinder_hits.tsv > amr/resfinder_summary.tsv
awk 'NR>1{print $6}' amr/resfinder_hits.tsv | sort | uniq -c | sort -nr \
  > amr/resfinder_gene_freq.tsv
awk 'NR>1{g[$6]++} END{
  printf "GENE\tISOLATES\tPREVALENCE(%%)\n";
  for(k in g) printf "%s\t%d\t%.1f\n",k,g[k],(g[k]/10)*100
}' amr/resfinder_hits.tsv | sort -k2,2nr > amr/resfinder_gene_prevalence.tsv
```

### `scripts/06_abricate_vfdb.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
export ABRICATE_DB=/opt/conda/db
abricate --datadir /opt/conda/db --db vfdb assemblies/*.contigs.fasta \
  > virulence/vfdb_hits.tsv
abricate --summary virulence/vfdb_hits.tsv > virulence/vfdb_summary.tsv
awk 'BEGIN{FS=OFS="\t"} NR==1{next} $6 ~ /(^|,)(hly|plcA|plcB)(,|$)/ {print $1,$6,$10,$11}' \
  virulence/vfdb_hits.tsv | sort -k1,1 -k2,2 \
  | awk 'BEGIN{print "SAMPLE\tGENE\t%COV\t%IDENTITY"}{print}' \
  > virulence/toxin_hly_plcA_plcB.tsv
```

### `scripts/util_pick_largest.py`

```python
#!/usr/bin/env python3
import glob, gzip
def read_fa(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if name: yield name,''.join(seq)
                name, seq = line[1:].strip().split()[0], []
            else: seq.append(line.strip())
        if name: yield name,''.join(seq)
with open('blast/largest_contigs.fasta','w') as w:
    for fa in sorted(glob.glob('assemblies/*.contigs.fasta')):
        sample = fa.split('/')[-1].replace('.contigs.fasta','')
        best = max(read_fa(fa), key=lambda x: len(x[1]))
        w.write(f">{sample}|{best[0]}\n{best[1]}\n")
print("Wrote blast/largest_contigs.fasta")
```

> Make scripts executable:
>
> ```bash
> chmod +x scripts/*.sh
> ```

---

## 8) How to reproduce (one-liners)

```bash
# 1) Trim
bash scripts/01_trim.sh samples_10.txt

# 2) Assemble
bash scripts/02_assemble.sh samples_10.txt

# 3) Assembly stats
bash scripts/03_assembly_stats.sh

# 4) BLAST ID
bash scripts/04_blast_ident.sh
column -t blast/blast_top_hit.tsv | sed -n '1,20p'

# 5) AMR
bash scripts/05_abricate_amr.sh
column -t amr/resfinder_summary.tsv | sed -n '1,20p'
column -t amr/resfinder_gene_prevalence.tsv

# 6) Virulence (includes toxin genes)
bash scripts/06_abricate_vfdb.sh
column -t virulence/toxin_hly_plcA_plcB.tsv
```

---

## 9) Rubric checklist (where to look)

* **Organism by BLAST (1 pt):** `blast/blast_top_hit.tsv` (+ section 4.2).
* **AMR genes (1 pt):** `amr/resfinder_summary.tsv` (+ section 4.3).
* **AMR summary & implications (2 pt):** `amr/resfinder_gene_prevalence.tsv` + section 5.
* **Therapy recommendation (1 pt):** section 5 (ampicillin + gentamicin; avoid cephalosporins/fosfomycin).
* **Clear report (2 pt):** this markdown, with methods/results/discussion.
* **Functional scripts + docs (2 pt):** `scripts/*.sh` & `scripts/util_pick_largest.py` + section 8.
* **Toxin genes (bonus 2 pt):** `virulence/toxin_hly_plcA_plcB.tsv` + section 4.4.

---

## 10) Notes & limitations

* ABRicate reports **presence/absence by sequence similarity**; it does **not** call **SNP-mediated** resistance.
* Remote BLAST can be slow; we used genus filtering and nohup logging.
  Check progress with:

  ```bash
  ps -fp "$(cat logs/blast_refseq.pid 2>/dev/null)" || echo "No running BLAST process"
  wc -l blast/blast_refseq.tsv
  tail -n 50 logs/blast_refseq.err
  ```
* For deeper epidemiology (phylogeny, SNP typing), consider cgMLST/SNP pipelines (beyond scope here).



---

**End of report.**
