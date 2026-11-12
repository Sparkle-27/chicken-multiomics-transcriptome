# A Genomic and Multi-tissue Transcriptomic Resource for Taihe Silky Chickens and Commercial Broilers

This repository provides analysis scripts and supporting information for a
multi-omics resource integrating whole-genome resequencing, short-read RNA
sequencing, and long-read transcriptome sequencing in chickens.

The resource was developed to support studies of genetic variation,
transcriptome annotation, isoform diversity, and gene expression across
metabolically relevant tissues. The main comparison focuses on slow-growing
Taihe Silky chickens and fast-growing commercial broilers. Publicly available
whole-genome sequencing data from additional chicken populations were also
included.

## 1. Study overview

The dataset comprises:

| Data type | Number | Description |
|---|---:|---|
| Whole-genome resequencing | 196 individuals | 100 newly sequenced individuals and 96 individuals from public datasets |
| Short-read RNA sequencing | 159 samples | Samples collected from four tissues |
| Long-read RNA sequencing | 24 samples | PacBio Iso-Seq data from four tissues |

The four tissues examined were:

| Abbreviation | Tissue |
|---|---|
| BM | Breast muscle |
| GM | Gastrocnemius muscle |
| LIV | Liver |
| VAF | Visceral abdominal fat |

The study design, sample selection, and data sources are described in the
associated manuscript and its supplementary tables.

## 2. Whole-genome sequencing data

### 2.1 Data generated in this study

Whole-genome resequencing data were newly generated for 100 chickens:

- 80 Taihe Silky chickens;
- 20 commercial broilers.

The Taihe Silky samples included individuals from the core cohort and an
additional cohort. Details of the cohorts and sample collection are provided
in the manuscript and Table 1.

### 2.2 Publicly available data

Previously published whole-genome resequencing data from 96 chickens were
included in the analysis:

| Population | Number of individuals | BioProject | Original publication |
|---|---:|---|---|
| Commercial broiler lines A and B | 39 | [PRJEB30270](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB30270) | [Qanbari et al., 2019](https://doi.org/10.1371/journal.pgen.1007989) |
| White Plymouth Rock | 57 | [PRJNA552722](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA552722) | [Guo et al., 2019](https://doi.org/10.3382/ps/pez411) |

These BioProjects contain samples beyond those included in the present study.
The specific samples used here are listed in the manuscript and associated
sample metadata.

## 3. Repository contents

The repository is organized into the following sections:

```text
.
├── README.md
├── scripts/
│   ├── wgs/
│   ├── short_read_rnaseq/
│   ├── long_read_isoseq/
│   ├── differential_expression/
│   └── figures/
├── config/
├── metadata/
└── results/
```

- `scripts/`: scripts used for data processing, analysis, and figure generation.
- `config/`: configuration files or parameter settings for the analyses.
- `metadata/`: sample information and links between sample identifiers and
  sequencing data.
- `results/`: selected processed results or instructions for accessing them.

Detailed descriptions of the scripts, input files, software versions, and
execution steps are provided in the corresponding directories.

## 4. Data availability

### 4.1 Newly generated sequencing data

The sequencing data generated in this study have been submitted to
[repository name].

- BioProject accession: **PRJCA065390**
- WGS data accession: **subSAM180981**
- Short-read RNA-seq data accession: **subSAM180499**
- PacBio Iso-Seq data accession: **subSAM182103**

### 4.2 Public sequencing data

The previously published whole-genome resequencing datasets are available from:

- Commercial broiler lines A and B:
  [BioProject PRJEB30270](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB30270).
- White Plymouth Rock chickens:
  [BioProject PRJNA552722](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA552722).

### 4.3 Processed data

Processed data and supporting files are available at:

- Sample metadata: **[add repository path or link]**
- Sequencing quality and mapping statistics: **[add repository path or link]**
- Transcript annotation, including GTF/GFF3 and FASTA files:
  **[add repository path or link]**
- Gene and transcript count matrices: **[add repository path or link]**
- Differential expression results: **[add repository path or link]**
- Other supporting files: **[add repository path or link]**

## 5. Reproducing the analyses

The analysis scripts are grouped by data type and analysis step. Before running
the scripts:

(1) Download the required sequencing data and reference files.  
(2) Install the software versions specified in the manuscript and script
   documentation.  
(3) Update input paths and analysis parameters in the relevant configuration
   files or scripts.  
(4) Run the scripts in the order described in each analysis directory.  

The reference genome assembly used in this study was **GRCg7b**, with
Ensembl annotation release 110. Additional reference files and analysis
parameters are specified in the manuscript and relevant scripts.

## 6. References for publicly available data

(1) Qanbari S, Rubin C-J, Maqbool K, et al. Genetics of adaptation in modern
   chicken. *PLoS Genetics*. 2019;15(4):e1007989.
   https://doi.org/10.1371/journal.pgen.1007989

(2) Guo Y, Lillie M, Zan Y, et al. A genomic inference of the White Plymouth Rock
   genealogy. *Poultry Science*. 2019;98(11):5272–5280.
   https://doi.org/10.3382/ps/pez411

## Citation

If you use this resource or the analysis scripts, please cite:

> **[Add the complete citation for the associated manuscript once available.]**

Please also cite the original publications when using the public datasets
listed above.

## Contact

For questions about the data or analyses, please contact:

**gonghuanfa@zju.edu.cn**  
**College of Animal Sciences, Zhejiang University, Hangzhou 310058, China**  
