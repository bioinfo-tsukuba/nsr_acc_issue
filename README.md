# Overview
This page serves as the supplement page for Hayakawa et al. in the NSR ACC issue.
A tutorial-style code is provided that aligns with the content of the paper using actual data.

The tutorial starts from loading scRNA-seq data in R and proceeds to a comparative analysis between wild-type and mutant samples.

The dataset used in this tutorial is from [Di Bella et al. 2021](https://www.nature.com/articles/s41586-021-03670-5).
In this tutorial, we use data from the cerebral cortex of the P1 stage of mouse development, comparing wild-type and Fezf2 mutant samples to identify genes that are specifically highly expressed in each condition.

# Tutorial

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

# How to run the tutorial on your own

## Get the code

1. Download this repositiory as a .zip file. [Download](https://github.com/bioinfo-tsukuba/nsr_acc_issue/archive/refs/heads/main.zip)
2. Uncompress `nsr_acc_issue-main.zip` to show the directory `nsr_acc_issue-main`.

## Preparation of scRNA-seq Data
### Overview of the dataset
The dataset consists of scRNA-seq data from the cerebral cortex at multiple developmental stages in wild-type and Fezf2 mutant mice, obtained from [Di Bella et al. 2021](https://www.nature.com/articles/s41586-021-03670-5).

### Data Sources

- [scRNA-seq Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153164) (GEO)
  1. In "Supplementary file" section, click "(http)" to download "GSE153164_RAW.tar".
  2. Uncompress the downloaded "GSE153164_RAW.tar".
  3. Put `GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5` and `GSM4635087_Fezf2KO_P1_filtered_feature_bc_matrix.h5` into `nsr_acc_issue-main/data/all`


## datasheet.tsv
A TSV file containing metadata for the scRNA-seq data, with relevant columns.

### Creating `datasheet.tsv`
This file includes columns specifying age, area, condition, and the path to the scRNA-seq data.
Then, put it into  `nsr_acc_issue-main/analysis/data_info`

#### Example `datasheet.tsv` [Link](https://github.com/bioinfo-tsukuba/nsr_acc_issue/blob/main/analysis/data_info/datasheet.tsv)
| age | area  | condition | path |
|-----|-------|-----------|-----------------------------------------|
| P1  | cortex | WT        | /home/rstudio/data/all/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5 |
| P1  | cortex | Fezf2KO   | /home/rstudio/data/all/GSM4635087_Fezf2KO_P1_filtered_feature_bc_matrix.h5 |

## Managing analysis environment using Docker
For this tutorial, we used an [Rstudio server Docker image](https://hub.docker.com/r/hway/rstudio_scrnaseq) that includes necessary libraries for scRNA-seq analysis.

# Additional information

## Cell Type Annotation Using Seurat FindTransferAnchors/DataTransfer

- In [04_celltype_annotation](https://bioinfo-tsukuba.github.io/nsr_acc_issue/04_celltype_annotation.html), by utilizing an annotated scRNA-seq dataset as a reference, cell types in the dataset of interest can be determined by mapping it to the reference data.
- Seurat provides functions (FindTransferAnchors / DataTransfer) specifically designed for cell type annotation using reference datasets. However, this method assumes that the cells in the reference dataset exhibit some degree of similarity to the cells in the dataset of interest. The annotation results depend on the cell types in the reference dataset. Additionally, if the gene names in the reference dataset differ significantly from those in the target dataset, this method cannot be applied. For this demonstration, the reference dataset used is the original annotated dataset from wild-type P1.
- The code to create the reference file is available from [mouse_p1_annotaion_prep.R](analysis/scripts/mouse_p1_annotaion_prep.R).
- The metadata used in `mouse_p1_annotaion_prep.R` can be prepared as follows:
  1. Go to the page on [Single Cell Portal - Broad Institute](https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex)
  2. Sign in (from the upper left button)
  3. Go to the page again
  4. Click "Download" to move to Download page
  5. Download "metaData_scDevSC.txt"
