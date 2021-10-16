
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/berkgurdamar/predatoR/blob/main/inst/extdata/predator_logo.png?raw=true" align="left" height=150/> predatoR: Mutation Impact Prediction Based on Network Properties

<!-- badges: start -->
<!-- badges: end -->

<br /> <br />

# Overview

`predatoR` is a tool for impact prediction of a mutation on a protein
structure by using network properties.

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
consist of *“PDB_ID”, “Chain”, “Position”, “Orig_AA”, “Mut_AA”*:

| PDB_ID | Chain | Position | Orig_AA | Mut_AA |
|:------:|:-----:|:--------:|:-------:|:------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |

Mutation impact prediction can be done via `predatoR()` function:

`predatoR()` function works on each PDB ID respectively. First,
downloads the PDB file, creates distance matrix and turns PDB structure
into a network. Calculates **Eigen Centrality Z-Score, Shortest Path
Z-score** and **Betweenness Z-Score** of input positions. Find the
**Gene Name** from the PDB file and if there are multiple genes
annotated for the input PDB file, asks user to choose the input gene,
gets **GnomAD Synonymous Z-Score, Non-Synonymous Z-Score, PLoF Score**,
and **Genic Intolerance Score**, gets **BLOSUM62 score** of the amino
acid change, finds the **KEGG Pathway Number** which contains the input
gene. Finally, make prediction based on an **Adaboost** model and
classifies the mutation as **Disease Causing** or **Silent**.

``` r
library(predatoR)
pred_res <- predatoR(input_df)
```

`predatoR()` function returns a data.frame which contains additional two
columns; **‘Prediction’** and **‘Probability’**. **‘Prediction’**
represents the result of the impact prediction and **‘Probability’**
represents the probability that the mutation classified as **Disease
Causing** or **Silent**.

| PDB_ID | Chain | Position | Orig_AA | Mut_AA | Prediction | Probability |
|:------:|:-----:|:--------:|:-------:|:------:|:----------:|:-----------:|
|  2DN2  |   B   |    1     |   VAL   |  ALA   |   Silent   |  0.6599836  |
|  2DN2  |   B   |    6     |   GLU   |  ALA   |   Silent   |  0.5899702  |
