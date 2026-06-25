# Acta2-LAD-remodeling: Multi-Omic Integration Pipeline

An advanced, end-to-end bioinformatic and statistical computing infrastructure integrating Bulk RNA-seq and CUT&Tag sequencing datasets. This repository is specifically engineered to investigate the dysregulation of **Lamin-Associated Domains (LADs)**, global **heterochromatin remodeling**, and aberrant transcriptional activation driven by the smooth muscle $\alpha$-actin (**Acta2**) **R179C** mutation in Vascular Smooth Muscle Cells (VSMCs).

---

## 1. 🔬 Scientific Background & Molecular Hypothesis

The smooth muscle $\alpha$-actin (*Acta2*) R179C missense mutation is clinically recognized as the primary genetic driver of a devastating, multi-systemic smooth muscle dysfunction syndrome. While alpha-actin is traditionally studied within the context of cytoplasmic microfilament contraction, recent cutting-edge structural and epigenetic profiling reveals a critical nuclear moonlighting role.

### The Mechanoeffect and LAD Collapse Hypothesis:
1. **Structural Disruption:** The R179C mutation alters actin polymerization dynamics, triggering a mechanical signaling failure that cascades into the nuclear envelope. This causes a structural collapse of nuclear lamins (specifically Lamin A/C and Lamin B1 structures).
2. **Epigenetic Unanchoring:** Under physiological conditions, transcriptionally repressed heterochromatin zones—known as **Lamin-Associated Domains (LADs)**—are securely anchored to the nuclear periphery. These domains are heavily enriched with repressive histone marks, primarily **H3K9me2** (Histone H3 Lysine 9 Dimethylation).
3. **Aberrant Activation:** The mutation-induced collapse of the lamin scaffold disrupts these anchoring points. Consequently, H3K9me2-marked heterochromatin undergoes a drastic genome-wide remodeling process. As these closed domains detach and untangle, previously silenced loci become pathologically exposed and active, driving an extensive, disease-specific transcriptional rewiring that governs severe vascular remodeling and phenotypic switching in VSMCs.

---

## 2. 📁 Comprehensive Repository Architecture

```text
Acta2-LAD-remodeling/
├── .gitattributes               # Git Large File Storage (LFS) tracking rules for heavy binary bigWig streams
├── README.md                    # Master documentation, scientific hypotheses, and deployment blueprint
├── CutTagR179C_WT.Rmd           # Main CUT&Tag script: Alignment parsing, duplication QC, and ChIPseeker mapping
├── PeakModifed.Rmd              # Window-based quantification pipeline: 100kb genomic tiling and KEGG pathways
├── cutTagFunction.R               # Production-grade R utility library (Overlap engines, species-independent ID mapping)
├── bulkRNAv3.html               # Compiled transcriptomic report validating original 3 reps with a 4th biological rep
└── h3k9me2_tracks_3rep/         # Specialized directory hosting high-resolution bigWig tracks formatted for browsers
    ├── R179C_Rep1_KWA-62279_S36_L001_bowtie2.mapped.bw
    ├── R179C_Rep2_KWA-62282_S39_L001_bowtie2.mapped.bw
    ├── R179C_Rep3_KWA-62285_S42_L001_bowtie2.mapped.bw
    ├── R179C_Group_Average.bw
    ├── wildType_Rep1_KWA-62277_S34_L001_bowtie2.mapped.bw
    ├── wildType_Rep2_KWA-62280_S37_L001_bowtie2.mapped.bw
    ├── wildType_Rep3_KWA-62283_S40_L001_bowtie2.mapped.bw
    ├── WT_Group_Average.bw
    ├── tgf_Rep1_KWA-62278_S35_L001_bowtie2.mapped.bw
    ├── tgf_Rep2_KWA-62281_S38_L001_bowtie2.mapped.bw
    ├── tgf_Rep3_KWA-62284_S41_L001_bowtie2.mapped.bw
    ├── R179C_vs_WT_Final_Log2Ratio.bw
    └── WT_vs_R179C_Final_Log2Ratio.bw
```

## 🛠️ Granular Pipeline Architectures & Code Implementation Details (bulkRNA.Rmd)

## 📊 Module A: Transcriptomic Cross-Validation & Pruning
To overcome standard sample size limitations and minimize false-discovery profiles, this workflow introduces a rigorous multi-tier verification layer:

Strict Expression Trimming: Preprocesses raw count matrices by implementing low-count and ultra-high-count filters. Outliers and background noise are trimmed out before DESeq2 dispersion modeling to maximize statistical sensitivity.

4th Biological Replicate Integration: Deploys a newly synchronized 4th biological replicate alongside the original 3-replicate core dataset to systematically isolate technical batch effects and evaluate cross-sample consistency.

## 🧬 Module B: Epigenomic Coordination Engine (CutTagR179C_WT.Rmd & cutTagFunction.R)
Manages downstream processing of H3K9me2 CUT&Tag distributions to decode heterochromatin patterns:

Consistency & Duplicate Metrics: Computes library duplication levels and tracks paired-end fragment size distribution profiles, ensuring high technical reproducibility across replicates.

ChIPseeker Annotation Matrix: Employs ChIPseeker::annotatePeak() coupled with TxDb.Mmusculus.UCSC.mm10.knownGene to assign peak coordinates to exact genomic traits. It is strictly configured to focus on Promoter regions defined within a flanking window of 5kb upstream to 1kb downstream of the Transcription Start Site (TSS -5kb to +1kb). It automatically handles cross-species identifier translations by mapping ENSEMBL vectors to NCBI Entrez IDs via clusterProfiler::bitr.

## 🧮 Module C: Gene-Centric 100kb Window Quantification (PeakModifed.Rmd)

Recognizing that LAD structures span broad, continuous megabase-scale genomic blocks rather than narrow, localized peaks, this pipeline bypasses standard peak-calling variants:Tiling Window Construction: Calculates the exact biological midpoint of every annotated gene body in the mouse genome. From each midpoint, it expands a standard 100kb genomic window ($\pm 50\text{ kb}$ upstream and downstream) using a safe boundary constraint:$$\text{Window Start} = \max(1, \text{Gene Center} - 50,000)$$$$\text{Window End} = \text{Gene Center} + 50,000$$Iterative Matrix Generation: Loops through compressed sample BED tracks (.bed.gz) and measures raw read overlaps within the 100kb boundaries via countOverlaps().Counts Per Million (CPM) Scaling: Normalizes raw matrices using calculated sequencing depths to permit comparative mathematical operations:$$\text{CPM} = \left( \frac{\text{Raw Counts}}{\text{Total Read Depth}} \right) \times 1,000,000$$Stringent Selection Metrics: Extracts mutation-specific candidates by demanding a multi-tier threshold check:Global Baseline Filter: Average global CPM must be $\ge 3$ to eliminate low-mappability regions.Enrichment Magnitude Check: Group Raw count average must be $\ge 150$.Directional Fold Change Constraint: Captures true R179C enrichment at $\log_2\text{FC} > 0.585$, and true WT enrichment at $\log_2\text{FC} < -0.585$.KEGG Network Profiling: Feeds verified Entrez ID blocks into clusterProfiler::enrichKEGG() under strict Benjamini-Hochberg false-discovery controls ($p \le 0.05$), rendering dotplots of pathways affected by LAD boundary adjustments.

## 🚀 Direct Load Configuration
Download the preconfigured tracking session document from this repository: R179CMutationGenomicBrwoser.json.

Launch the official WashU Epigenome Browser Gateway.
