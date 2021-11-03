
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

<!-- <br /> -->
<!-- <br /> -->

# Overview

`predatoR` is a tool for impact prediction of a mutation on a protein
structure by using network properties.

<img src="https://github.com/berkgurdamar/predatoR/blob/main/vignettes/predatoR_workflow.png?raw=true" style="max-width:100%;" />

`predatoR()` function works on each PDB ID respectively. For each PDB
ID;

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
**‘Mut_AA’** and **‘Gene_Name’** as an optional input:

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |    HBB    |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |    HBB    |

Mutation impact prediction can be done via `predatoR()` function:

`predatoR()` function works on each PDB ID respectively. First,
downloads the PDB file, creates distance matrix between every atom and
turns PDB structure into a network model. Calculates **Degree Centrality
Z-Score, Eigen Centrality Z-Score, Shortest Path Z-score, Betweenness
Z-Score, Clique Z-Score** and **PageRank Z-Score** of input positions.
If **‘Gene_Name’** is not provided in the input, `predatoR()` gets the
related gene names from **‘Ensembl’** and if there are multiple genes
annotated for the same PDB ID, gets the gene that has maximum
**Synonymous Z-Score, Non-Synonymous Z-Score**, and **PLoF Score**, and
**Genic Intolerance Score**, gets **BLOSUM62 score** of the mutation,
finds the **KEGG Pathway Number** which contains the input gene.
Finally, make prediction based on a pre-computed **Adaboost** model and
classifies the mutation as **Disease Causing** or **Silent**.

``` r
library(predatoR)

input_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

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
