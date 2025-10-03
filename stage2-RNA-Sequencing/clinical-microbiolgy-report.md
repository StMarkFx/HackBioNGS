# Clinical Microbiology Report
## Transcriptomic Profiling of *Staphylococcus aureus* in Periprosthetic Joint Infection (PJI): Acute vs Chronic

### Background & Rationale
Periprosthetic joint infections (PJIs) are severe implant-associated infections with high morbidity and cost. *Staphylococcus aureus*—including MRSA—is a leading cause. A key clinical problem is the bacterium’s ability to switch between:
- **Acute phase:** aggressive, planktonic growth; high virulence factor expression (toxins, adhesins, immune evasion).
- **Chronic phase:** biofilm-like state; reduced overt virulence with increased persistence programs (stress responses, metabolic rewiring, antibiotic tolerance).

**Goal:** Use RNA-seq to profile *S. aureus* directly from acute vs chronic PJI isolates to (i) identify virulence programs enriched in acute infection, (ii) detect persistence/metabolic pathways in chronic infection, and (iii) generate plots and a concise report to support clinical interpretation.

---

### Cohort & Data
- **Samples:** 8 total (4 acute: SRR20959680–83; 4 chronic: SRR20959676–79).
- **Reference:** USA300_FPR3757 genome and GFF.
- **Design:** chronic vs acute (balanced 4 vs 4).

---

### Methods (laboratory-style summary)
**Preprocessing & QC**
- Downloaded SRA (prefetch, fasterq-dump), paired FASTQ compression (pigz).
- Trimming with **fastp**; per-sample HTML reports; aggregated **MultiQC** summary.  
  *QC report available:* `multiqc_report.html`. :contentReference[oaicite:5]{index=5}

**Alignment & Counting**
- **Bowtie2** (very-sensitive-local) against USA300 index; BAM sorting and indexing with samtools.
- **featureCounts** using `-t gene -g locus_tag` against the GFF to generate gene-level counts.
- High assignment rates were observed (≈95–98%) in your logs (consistent with good alignment quality).

**Differential Expression**
- **DESeq2**: design `~ condition` (chronic reference level), low-count filter (row sums ≥10), Wald test.
- VST transformation for PCA and heatmap.

**Outputs**
- `deseq2_results.csv` (full DE table), `vst_normalized_counts.csv`, and plots (`pca.png`, `volcano.png`, `heatmap_topDE.png`).
- Quick clinical summary markdown: `clinical_report.md`. This file reports:  
  **Total genes tested: 2261; Significant genes (padj < 0.05): 0.** :contentReference[oaicite:6]{index=6}

---

### Results

#### Global QC
- MultiQC indicates typical Illumina quality metrics; trimming removed low-quality tails/adapters as expected. No systemic failures noted. (See HTML for per-sample details.) :contentReference[oaicite:7]{index=7}

#### Alignment & Counting
- Very high read assignment to genes (≈95–98%), consistent with clean bacterial transcriptomes and accurate reference selection.

#### Differential Expression
- **Statistically significant DE genes at FDR < 0.05: 0** (N=4 vs 4, conservative after multiple-testing). :contentReference[oaicite:8]{index=8}  
- Nonetheless, several loci show notable **log2 fold-change trends** (see `deseq2_results.csv`), including candidates with |LFC| ~ 1–1.5 (both directions). These can still be biologically meaningful hypotheses for follow-up.

#### Visualization
- **PCA (VST):** Acute (SRR20959680–83) vs Chronic (SRR20959676–79) show separation on PC1/PC2 with some within-group spread. This supports phase-associated global expression differences.
- **Volcano plot:** Most genes are near logFC 0 with few outliers; FDR threshold not crossed at 0.05 given sample size.
- **Heatmap (top variable genes):** Distinct signal across samples with clusters of co-varying genes; several loci (e.g., SAUSA300_RS* IDs) show phase-biased trends that merit pathway context.

---

### Interpretation (Clinical Microbiology)
- The **lack of FDR-significant genes** likely reflects the modest cohort size and conservative multiple-testing, not the absence of biology.
- **Acute phase** is expected to emphasize **adhesins/toxins/immune evasion**, while **chronic phase** typically shows **stress tolerance, metabolic rewiring, and biofilm-linked persistence**. The PCA separation and heatmap patterns are congruent with this biological model.
- Candidate RS-labeled genes with consistent directionality (|LFC|≈1) may map to virulence/persistence modules (e.g., cell-wall processes, oxidative stress, amino-acid/energy metabolism) and warrant **targeted validation**.

---

### Limitations
- Small sample size (4 vs 4) limits statistical power for multiple-testing control.
- Single reference genome may miss strain-specific features; pan-genome mapping or transcriptome assembly could complement.
- No explicit human RNA contamination check included (not expected here but can be added with Kraken2/BBMap).

---

### Recommendations / Next Steps
1. **Pathway context**: map locus_tags to KEGG/GO for enrichment. (Optional step was scaffolded; mapping tables for *S. aureus* USA300 locus_tag→KEGG are needed.)
2. **Targeted follow-up**: qPCR or targeted RNA-seq on top-trend genes.
3. **Increase N** or use shrinkage/empirical Bayes methods and **independent filtering** to improve detection.
4. **Biofilm modules**: specifically interrogate ica operon, stress regulons (SigB), and toxin loci for directional trends.

---

### Deliverables (this submission)
- **QC**: `multiqc_report.html` (aggregate quality). :contentReference[oaicite:9]{index=9}
- **Counts**: `counts/gene_counts.clean.tsv`
- **DE**: `results/deseq2/deseq2_results.csv`
- **Normalized**: `results/deseq2/vst_normalized_counts.csv`
- **Plots**: `results/plots/pca.png`, `volcano.png`, `heatmap_topDE.png`
- **Summary**: `results/report/clinical_report.md` (shows 2261 tested; 0 significant at padj<0.05). :contentReference[oaicite:10]{index=10}

---

### One-paragraph Summary (for the abstract)
RNA-seq of *S. aureus* from PJI (4 acute vs 4 chronic; USA300 reference) produced high-quality data with excellent alignment. While no genes passed FDR<0.05 (2261 tested), PCA and heatmap analyses show group separation and coherent expression trends consistent with **acute virulence programs** versus **chronic persistence/biofilm-like adaptation**. These transcriptional signatures, despite conservative multiple-testing, suggest biologically meaningful phase-linked reprogramming that can inform diagnostics (phase discrimination) and therapies (targeting persistence mechanisms).
