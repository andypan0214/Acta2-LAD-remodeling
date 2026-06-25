# Acta2-LAD-remodeling

Multi-omic data integration pipeline (RNA-seq & CUT&Tag) investigating epigenetic mechanisms in *Acta2* R179C/+ smooth muscle cells.

---

## 📊 Project Overview
The smooth muscle $\alpha$-actin (**Acta2**) **R179C** mutation is clinically associated with severe smooth muscle dysfunction and multi-systemic smooth muscle dysfunction syndrome. This repository hosts the centralized analysis infrastructure to explore how this specific mutation alters the structural organization of **Lamin-Associated Domains (LADs)** and triggers genome-wide heterochromatin remodeling.

## 📁 Repository Infrastructure
```text
Acta2-LAD-remodeling/
├── .gitattributes               # Git LFS configuration managing bigWig tracking
├── README.md                    # Project documentation & visualization configuration
├── CutTagR179C_WT.Rmd           # Master pipeline: QC metrics, ChIPseeker, and RNA-seq cross-referencing
└── h3k9me2_tracks_3rep/         # High-resolution bigWig alignment files (3 biological replicates)

