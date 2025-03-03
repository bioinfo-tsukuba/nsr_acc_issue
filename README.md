# nsr_acc_issue

# Overview
This page serves as the supplement page for Hayakawa et al. in the NSR ACC issue.
A tutorial-style code is provided that aligns with the content of the paper using actual data.

The tutorial starts from loading scRNA-seq data in R and proceeds to a comparative analysis between wild-type and mutant samples.

The dataset used in this tutorial is from [Di Bella et al. 2021](https://www.nature.com/articles/s41586-021-03670-5).
In this tutorial, we use data from the cerebral cortex of the P1 stage of mouse development, comparing wild-type and Fezf2 mutant samples to identify genes that are specifically highly expressed in each condition.

The analysis consists of the following five steps, each linked to a specific code implementation:

- [01_load](https://bioinfo-tsukuba.github.io/nsr_acc_issue/01_load.html)
  - Load scRNA-seq data in R and prepare it for analysis with Seurat.
- [02_quality_control](https://bioinfo-tsukuba.github.io/nsr_acc_issue/02_qaulity_control.html)
  - Perform quality control on scRNA-seq data and remove low-quality cells.
- [03_data_preprocessing_and_integration](https://bioinfo-tsukuba.github.io/nsr_acc_issue/03_data_preprocessing_and_integration.html)
  - Execute the basic workflow for scRNA-seq data analysis.
- [04_celltype_annotation](https://bioinfo-tsukuba.github.io/nsr_acc_issue/04_celltype_annotation.html)
  - Annotate cell types for each cell. This section introduces three different methods for cell type annotation.
- [05_comparative_analysis](https://bioinfo-tsukuba.github.io/nsr_acc_issue/05_comparative_analysis.html)
  - Perform a comparative analysis between wild-type and Fezf2 mutants to identify genes that are specifically highly expressed in each condition.

# Additional Information
## Preparation of scRNA-seq Data
### Overview of the Dataset
The dataset consists of scRNA-seq data from the cerebral cortex at multiple developmental stages in wild-type and Fezf2 mutant mice, obtained from [Di Bella et al. 2021](https://www.nature.com/articles/s41586-021-03670-5).

### Data Sources
- [scRNA-seq Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153164) (GEO)
- [Metadata](https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex)

## datasheet.tsv
A TSV file containing metadata for the scRNA-seq data, with relevant columns.

### Creating `datasheet.tsv`
This file includes columns specifying age, area, condition, and the path to the scRNA-seq data.

#### Example `datasheet.tsv` [Link](https://github.com/bioinfo-tsukuba/nsr_acc_issue/blob/main/analysis/data_info/datasheet.tsv)
| age | area  | condition | path |
|-----|-------|-----------|-----------------------------------------|
| P1  | cortex | WT        | /home/rstudio/data/all/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5 |
| P1  | cortex | Fezf2KO   | /home/rstudio/data/all/GSM4635087_Fezf2KO_P1_filtered_feature_bc_matrix.h5 |

## Managing Analysis Environments Using Docker
For this tutorial, we used an [Rstudio server Docker image](https://hub.docker.com/r/hway/rstudio_scrnaseq) that includes necessary libraries for scRNA-seq analysis.
