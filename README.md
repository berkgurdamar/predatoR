
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

<!-- <br /> -->
<!-- <br /> -->

# Overview

`predatoR` is a tool for mutation impact prediction based on network
properties.

<img src="https://github.com/berkgurdamar/predatoR/blob/main/vignettes/predatoR_workflow.png?raw=true" style="max-width:100%;" />

`predatoR()` function is the wrapper function of `predatoR` package.

`predatoR()` works on each PDB respectively. For each PDB;

-   Download/Read PDB file
-   Calculate distances between every atom
-   Create network from PDB
-   Calculate Degree Centrality Z-Score of each atom
-   Calculate Eigen Centrality Z-Score of each atom
-   Calculate Shortest Path Z-Score of each atom
-   Calculate Betweenness Z-Score of each atom
-   Calculate Clique Z-Score of each atom
-   Calculate PageRank Z-Score of each atom
-   Gets gnomAD Synonymous Z-Score, Non-Synonymous Z-Score, and PLoF
    Score
-   Gets BLOSUM62 score of the mutation
-   Finds the number of KEGG Pathways which contains the input gene
-   Gets [Genic Intolerance](http://genic-intolerance.org/) Score
-   Finds the number of GO terms associated with the input gene
-   Finds the number of diseases associated with the input gene from
    [DisGeNET](https://www.disgenet.org/)
-   Gets the Gene Essentiality Score from Online Gene Essentiality
    ([OGEE](https://v3.ogee.info/#/home)) Database
-   Calculates the median gene expression value of the input gene from
    [GTEx](https://gtexportal.org/home/)
-   Calculates 6 different features from Accessible Surface Area and
    Hydrophobicity Scale of reference and mutant amino acids
-   Make prediction

# Installation

You can install predatoR via
[devtools](https://www.r-project.org/nosvn/pandoc/devtools.html):

``` r
library(devtools)
install_github("berkgurdamar/predatoR")
```

# Usage

Mutation impact prediction can be done via `predatoR()` function:

`predatoR()` uses data.frame structures as an input. data.frame should
consist of **‘PDB_ID’**, **‘Chain’**, **‘Position’**, **‘Orig_AA’**,
**‘Mut_AA’** and **‘Gene_Name’** (optional):

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|
|  4RFZ  |   A   |   414    |   GLY   |  ARG   |    BTK    |
|  1Z2M  |   A   |    21    |   SER   |  ASN   |   ISG15   |

`predatoR()` can work with input which has partially included gene
names.

``` r
library(predatoR)

# Gene name included
input_df <- as.data.frame(rbind(c("4RFZ", "A", 414, "GLY", "ARG", "BTK"),
                                c("1Z2M", "A", 21,  "SER", "ASN", "ISG15")))

pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE)

# Gene name not included
input_df <- as.data.frame(rbind(c("4RFZ", "A", 414, "GLY", "ARG"),
                                c("1Z2M", "A", 21,  "SER", "ASN")))

pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = FALSE)

# Partially included gene names
input_df <- as.data.frame(rbind(c("4RFZ", "A", 414, "GLY", "ARG", "BTK"),
                                c("1Z2M", "A", 21,  "SER", "ASN", "")))

pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE)
```

`predatoR()` function returns a data.frame which contains additional two
columns; **‘Prediction’** and **‘Probability’**. **‘Prediction’**
represents the result of the impact prediction and **‘Probability’**
represents the probability that the mutation classified as
**Pathogenic** or **Neutral**.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name | Prediction | Probability |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|:----------:|:-----------:|
|  4RFZ  |   A   |   414    |   GLY   |  ARG   |    BTK    | Pathogenic |  0.7666323  |
|  1Z2M  |   A   |    21    |   SER   |  ASN   |   ISG15   |  Nautral   |  0.7704317  |

## Utility Functions

The wrapper function `predatoR()` uses the utility functions below;

-   `read_PDB()`
-   `PDB2connections()`
-   `degree_score()`
-   `eigen_centrality_score()`
-   `shorteset_path_score()`
-   `betweenness_score()`
-   `clique_score()`
-   `pagerank_score()`
-   `gnomad_scores()`
-   `BLOSUM62_score()`
-   `KEGG_pathway_number()`
-   `genic_intolerance()`
-   `GO_terms()`
-   `DisGeNET()`
-   `gene_essentiality()`
-   `GTEx()`
-   `amino_acid_features()`
-   `impact_prediction()`

Utility functions can be used alone, for more detail please see
vignette.
