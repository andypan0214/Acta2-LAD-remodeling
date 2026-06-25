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
├── cutTag_utils.R               # Production-grade R utility library (Overlap engines, species-independent ID mapping)
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

