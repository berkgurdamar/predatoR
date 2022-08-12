
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/berkgurdamar/predatoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/berkgurdamar/predatoR/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/berkgurdamar/predatoR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/berkgurdamar/predatoR?branch=main)
<!-- badges: end -->

<!-- <br /> -->
<!-- <br /> -->
<!-- # NOTE -->
<!-- `R-CMD-check` and `Codecov` gave an error due to the Git LFS. Local tests are fine you can download the `predatoR` as described below. -->

# Overview

`predatoR` is a tool for mutation impact prediction based on network
properties.

<img src="https://github.com/berkgurdamar/predatoR/blob/main/vignettes/predatoR_workflow.png?raw=true" style="max-width:100%;" />

`predatoR()` function is the wrapper function of `predatoR` package.

`predatoR()` works on each PDB respectively. For each PDB;

-   Downloads/Reads PDB file
-   Calculates distance between each atom in the structure
-   Creates a network from PDB
-   Calculates Degree Centrality Z-Score of each atom
-   Calculates Eigen Centrality Z-Score of each atom
-   Calculates Shortest Path Z-Score of each atom
-   Calculates Betweenness Z-Score of each atom
-   Calculates Clique Z-Score of each atom
-   Calculates PageRank Z-Score of each atom
-   Gets [gnomAD](https://gnomad.broadinstitute.org/) Synonymous
    Z-Score, Missense Z-Score, and pLI Score
-   Gets BLOSUM62 score of the mutation
-   Finds the number of [KEGG Pathways](https://www.genome.jp/kegg/)
    which contains the input gene
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
-   Makes prediction

# Installation

You can install the predatoR via
[devtools](https://www.r-project.org/nosvn/pandoc/devtools.html):

``` r
library(devtools)
install_github("berkgurdamar/predatoR")
```

# Usage

Mutation impact prediction can be done via `predatoR()` function:

`predatoR()` uses data.frame structures as an input. data.frame should
consist of **‘PDB_ID’**, **‘Chain’**, **‘Position’**, **‘Orig_AA’**,
**‘Mut_AA’** and **‘Gene_Name’** (optional). Predictions can be made by
using 2 different models, 7 Angstrom (Å)-all atoms model and 7Å-carbon
alpha (Cα) atoms only model.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|
|  3SQJ  |   A   |   196    |   GLN   |  LEU   |    ALB    |
|  3SQJ  |   A   |   396    |   GLU   |  LYS   |    ALB    |

`predatoR()` can work with input which has partially included gene
names.

``` r
library(predatoR)
# Gene name included
input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU", "ALB"),
                               c("3SQJ", "A", 396, "GLU", "LYS", "ALB")))
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE)
# Gene name not included
input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU"),
                               c("3SQJ", "A", 396, "GLU", "LYS")))
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = FALSE)
# Partially included gene names
input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU", "ALB"),
                               c("3SQJ", "A", 396, "GLU", "LYS", "")))
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE)
```

`predatoR()` function returns a data.frame which contains additional two
columns; **‘Prediction’** and **‘Probability’**. **‘Prediction’**
represents the result of the impact prediction and **‘Probability’**
represents the probability that the mutation classified as
**Pathogenic** or **Neutral**.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name | Prediction | Probability |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|:----------:|:-----------:|
|  3SQJ  |   A   |   196    |   GLN   |  LEU   |    ALB    |  Neutral   |  0.6521832  |
|  3SQJ  |   A   |   396    |   GLU   |  LYS   |    ALB    |  Neutral   |  0.6009792  |

## Exploratory Usage

Network properties can be calculated by using different distance
cutoffs. In this approach, `predatoR()` does not make any prediction
about the mutation, but returns a data frame contains all 24 features
annotated to the dataset. Both network formalisation approaches also can
be used.

``` r
# networks build using all atoms and 7.6Å cutoff
prediction_result <- predatoR(input_df, distance_cutoff = 7.6, network_approach = "all") 

# networks build using only Cα atoms and 8Å cutoff
prediction_result <- predatoR(input_df, distance_cutoff = 8, network_approach = "ca") 
```

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
