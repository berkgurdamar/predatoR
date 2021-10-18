
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->
<!-- badges: end -->

<br /> <br />

# Overview

`predatoR` is a tool for impact prediction of a mutation on a protein
structure by using network properties.

<img src="https://github.com/berkgurdamar/predatoR/blob/main/vignettes/predatoR_workflow.png?raw=true" style="max-width:100%;" />

`predatoR()` function works on each PDB ID respectively. For each PDB
ID;

-   Download/Read PDB file
-   Calculate distances between every atom
-   Create network from PDB
-   Calculate Eigen Centrality Z-Score of each Carbon-α atom
-   Calculate Shortest Path Z-Score of each Carbon-α atom
-   Calculate Betweenness Z-Score of each Carbon-α atom
-   Gets GnomAD Synonymous Z-Score, Non-Synonymous Z-Score, and PLoF
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
turns PDB structure into a network model. Calculates **Eigen Centrality
Z-Score, Shortest Path Z-score** and **Betweenness Z-Score** of input
positions. If **‘Gene_Name’** is not provided in the input, `predatoR()`
gets the related gene names from **‘Ensembl’** and if there are multiple
genes annotated for the PDB ID, asks user to choose input
**‘Gene_Name’**, gets **GnomAD Synonymous Z-Score, Non-Synonymous
Z-Score, PLoF Score**, and **Genic Intolerance Score**, gets **BLOSUM62
score** of the mutation, finds the **KEGG Pathway Number** which
contains the input gene. Finally, make prediction based on a
pre-computed **Adaboost** model and classifies the mutation as **Disease
Causing** or **Silent**.

``` r
library(predatoR)
pred_res <- predatoR(input_df)
```

`predatoR()` function returns a data.frame which contains additional two
columns; **‘Prediction’** and **‘Probability’**. **‘Prediction’**
represents the result of the impact prediction and **‘Probability’**
represents the probability that the mutation classified as **Disease
Causing** or **Silent**.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Gene_Name | Prediction | Probability |
|:------:|:-----:|:--------:|:-------:|:------:|:---------:|:----------:|:-----------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |    HBB    |   Silent   |  0.6205009  |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |    HBB    |   Silent   |  0.6286857  |
