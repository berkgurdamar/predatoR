
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

<!-- <br /> -->
<!-- <br /> -->

# Overview

`predatoR` is a tool for mutation impact prediction using network
properties.

<img src="https://github.com/berkgurdamar/predatoR/blob/main/vignettes/predatoR_workflow.png?raw=true" style="max-width:100%;" />

`predatoR()` function is the wrapper function of `predatoR` package.

`predatoR()` works on each PDB ID respectively. For each PDB ID;

-   Download/Read PDB file
-   Calculate distances between every atom
-   Create network from PDB
-   Calculate Degree Centrality Z-Score of each Carbon-α atom
-   Calculate Eigen Centrality Z-Score of each Carbon-α atom
-   Calculate Shortest Path Z-Score of each Carbon-α atom
-   Calculate Betweenness Z-Score of each Carbon-α atom
-   Calculate Clique Z-Score of each Carbon-α atom
-   Calculate PageRank Z-Score of each Carbon-α atom
-   Gets gnomAD Synonymous Z-Score, Non-Synonymous Z-Score, and PLoF
    Score
-   Gets Genic Intolerance Score
-   Gets BLOSUM62 score of the mutation
-   Finds the number of KEGG Pathways which contains the input gene
-   Make prediction

# Installation

You can install the released version of predatoR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("predatoR")
```

or install via
[devtools](https://www.r-project.org/nosvn/pandoc/devtools.html):

``` r
library(devtools)
install_github("berkgurdamar/predatoR")
```

# Usage

`predatoR` uses data.frame structures as an input. data.frame should
consist of **‘PDB_ID’**, **‘Chain’**, **‘Position’**, **‘Orig_AA’**,
**‘Mut_AA’** and **‘Gene_Name’**(optional):

Mutation impact prediction can be done via `predatoR()` function:

``` r
input_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                c("2DN2", "B", 6, "GLU", "ALA", "HBB")))
```

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |    HBB    |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |    HBB    |

``` r
library(predatoR)

pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE)
```

`predatoR()` function returns a data.frame which contains additional two
columns; **‘Prediction’** and **‘Probability’**. **‘Prediction’**
represents the result of the impact prediction and **‘Probability’**
represents the probability that the mutation classified as **Disease
Causing** or **Silent**.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name | Prediction | Probability |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|:----------:|:-----------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |    HBB    |   Silent   |  0.6123515  |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |    HBB    |   Silent   |  0.5398617  |

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
-   `impact_prediction()`

Utility functions can be used alone, for more detail please see
vignette.
