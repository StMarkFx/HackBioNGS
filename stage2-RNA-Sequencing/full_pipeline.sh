1) Reproducible pipeline script (save as scripts/00_full_pipeline.sh)

Copy the block below into scripts/00_full_pipeline.sh, then:

chmod +x scripts/00_full_pipeline.sh
tmux new -s pji   # optional, so long steps keep running
./scripts/00_full_pipeline.sh

#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# PJI RNA-seq: Acute vs Chronic S. aureus (PRJNA867318)
# End-to-end, commented pipeline reproducing the work you ran.
#
# NOTE ABOUT LONG-RUNNING STEPS:
#  - Marked with [LONG]. Run inside tmux (Ctrl+b then d to detach).
#  - Monitor with: tail -f logs/<something>.log
# ============================================================

# --------- 0) Project layout ----------
# (Already created in your run; safe to re-run)
mkdir -p ~/pji-transcriptome
cd       ~/pji-transcriptome

mkdir -p envs config \
  refs/{genomes,annotations} \
  data/{sra,fastq} \
  qc/{fastqc,multiqc} \
  trim align/{bam,logs} counts \
  results/{deseq2,plots,enrichment,report} \
  scripts logs

# --------- 1) Conda env (tools + R stack) ----------
# You installed an env named pji-rnaseq with sra-tools, fastp, fastqc, bowtie2,
# samtools, subread, pigz, yq, R (>=4.3), tidyverse, DESeq2, pheatmap, ggrepel, etc.
# If recreating on a fresh machine, put your env YAML at envs/env.yaml and do:
# mamba env create -f envs/env.yaml || conda env create -f envs/env.yaml
# echo "conda activate pji-rnaseq" > activate.sh
source ./activate.sh   # activates the env you already created

# --------- 2) SRA download ----------
# Accessions (4 chronic, 4 acute)
cat > data/sra/acc.txt <<'ACC'
SRR20959676
SRR20959677
SRR20959678
SRR20959679
SRR20959680
SRR20959681
SRR20959682
SRR20959683
ACC

# Keep SRA cache local so it’s easy to clean up if needed
vdb-config --prefetch-to-cwd || true

# [LONG] Download .sra
# Monitor: tail -f logs/prefetch.log
cat data/sra/acc.txt | parallel -j 4 "prefetch {} -O data/sra" |& tee logs/prefetch.log

# [LONG] Convert to paired FASTQ, then gzip
# Monitor: tail -f logs/fasterq.log
cat data/sra/acc.txt | parallel -j 2 "
  fasterq-dump --split-files --threads 8 -O data/fastq data/sra/{}
  pigz -p 8 data/fastq/{}_*.fastq
" |& tee -a logs/fasterq.log

# --------- 3) QC + trimming ----------
# FastQC (raw)
fastqc -t 8 -o qc/fastqc data/fastq/*.fastq.gz |& tee logs/fastqc_raw.log

# Trimming (adapter/quality) with fastp [fast, but many files]
# Creates per-sample HTML/JSON reports in qc/fastqc
for fq1 in data/fastq/*_1.fastq.gz; do
  fq2=${fq1/_1.fastq.gz/_2.fastq.gz}
  sample=$(basename "$fq1" | sed 's/_1\.fastq\.gz//')
  fastp -i "$fq1" -I "$fq2" \
        -o "trim/${sample}_1.trim.fastq.gz" -O "trim/${sample}_2.trim.fastq.gz" \
        --thread 8 --detect_adapter_for_pe \
        --html "qc/fastqc/${sample}_fastp.html" \
        --json "qc/fastqc/${sample}_fastp.json" \
        |& tee -a logs/fastp.log
done

# FastQC (trimmed) and MultiQC summary (open this to inspect quality)
fastqc -t 8 -o qc/fastqc trim/*.trim.fastq.gz |& tee logs/fastqc_trimmed.log
multiqc qc -o qc/multiqc |& tee logs/multiqc.log
# The file qc/multiqc/multiqc_report.html aggregates all QC. :contentReference[oaicite:1]{index=1}

# --------- 4) Reference download + index (USA300_FPR3757) ----------
# Genome + GFF (if not already present)
# (Your run already placed these in refs/, and Bowtie2 index exists.)
# Example:
# curl -L -o refs/genomes/genomic.fna.gz <NCBI-FASTA-URL>
# curl -L -o refs/annotations/genomic.gff.gz <NCBI-GFF-URL>
# gunzip -c refs/genomes/genomic.fna.gz > refs/genomes/genomic.fna
# gunzip -c refs/annotations/genomic.gff.gz > refs/annotations/genomic.gff

# Build Bowtie2 index (already done; safe to skip if 6 files exist)
if [ ! -f refs/genomes/saureus_bt2.1.bt2 ]; then
  bowtie2-build refs/genomes/genomic.fna refs/genomes/saureus_bt2 |& tee logs/bt2_build.log
fi

# --------- 5) Alignment (Bowtie2) + sorting (samtools) ----------
# [LONG] Align each pair; write sorted BAM + index
for fq1 in trim/*_1.trim.fastq.gz; do
  fq2=${fq1/_1.trim.fastq.gz/_2.trim.fastq.gz}
  sample=$(basename "$fq1" | sed 's/_1\.trim\.fastq\.gz//')
  bowtie2 --very-sensitive-local -x refs/genomes/saureus_bt2 -1 "$fq1" -2 "$fq2" -p 8 \
    2> "align/logs/${sample}.bt2.log" \
  | samtools view -bS - \
  | samtools sort -@ 4 -o "align/bam/${sample}.sorted.bam"
  samtools index "align/bam/${sample}.sorted.bam"
done

# --------- 6) Gene counting (featureCounts) ----------
# -t gene -g locus_tag uses the gene features and locus tags from USA300 GFF
featureCounts -T 8 -p -B -C \
  -t gene -g locus_tag \
  -a refs/annotations/genomic.gff \
  -o counts/gene_counts.txt align/bam/*.sorted.bam |& tee logs/featureCounts.log

# Clean sample headers (drop “align/bam/” and “.sorted.bam”)
awk 'BEGIN{FS=OFS="\t"} NR==2{for(i=7;i<=NF;i++) gsub(/.*\//,"",$i); for(i=7;i<=NF;i++) gsub(/\.sorted\.bam/,"",$i); print; next} NR>2{print}' \
  counts/gene_counts.txt > counts/gene_counts.clean.tsv

# --------- 7) DESeq2 (R) ----------
# Uses your existing scripts/06_dge.R which:
#  - reads counts/gene_counts.clean.tsv
#  - builds chronic vs acute design
#  - runs DESeq2, saves results table and VST matrix
Rscript scripts/06_dge.R |& tee logs/deseq2.log

# Results produced:
#   results/deseq2/deseq2_results.csv     (full DE table)
#   results/deseq2/vst_normalized_counts.csv
#   results/deseq2/{dds.rda,res_df.rda}   (R objects)
#   results/plots/{pca.png,volcano.png,heatmap_topDE.png}

# --------- 8) Optional enrichment (when mapping is available) ----------
# Rscript scripts/07_enrichment.R |& tee logs/enrichment.log  # (optional)

# --------- 9) Plots (already produced by your scripts) ----------
# Rscript scripts/08_plots.R       |& tee logs/plots.log
# Rscript scripts/08_plots_fix.R   |& tee logs/plots_fix.log

# --------- 10) Report generation ----------
# You created a markdown summary:
#   results/report/clinical_report.md  (quick clinical-style summary) :contentReference[oaicite:2]{index=2}
# Optionally render a full HTML/PDF with RMarkdown if pandoc/tinytex are installed.

echo "Pipeline complete.
Key files:
  qc/multiqc/multiqc_report.html
  counts/gene_counts.clean.tsv
  results/deseq2/deseq2_results.csv
  results/deseq2/vst_normalized_counts.csv
  results/plots/pca.png
  results/plots/volcano.png
  results/plots/heatmap_topDE.png
  results/report/clinical_report.md
"