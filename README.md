
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PredImption: Mutation Impact Predition Based on Network Properties

<!-- badges: start -->
<!-- badges: end -->

# Overview

‘PredImption’ is a tool for impact prediction of a mutation on a protein
structure by using network properties.

# Installation

You can install the released version of PredImption from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PredImption")
```

or install via ‘devtools’:

``` r
library(devtools)
devtools::install_github("berkgurdamar/PredImption")
```

# Usage

‘PredImption’ uses data.frame structures as an input. data.frame should
consist “PDB_ID”, “Chain”, “Position”, “Orig_AA”, “Mut_AA”:

|  V1  | V2  | V3  | V4  | V5  |
|:----:|:---:|:---:|:---:|:---:|
| 2DN2 |  B  |  1  | VAL | ALA |
| 2DN2 |  B  |  6  | GLU | ALA |
| 4ONL |  A  | 36  | GLU | GLY |
| 4ONL |  A  | 78  | PRO | GLN |

Mutation impact prediction can be done via `PredImption()` function:

``` r
library(PredImption)
pred_res <- PredImption(input_df)
```

`PredImption()` function works on each PDB ID respectively. First, gets
the PDB structure, create distance matrix and turns PDB structure into a
network. Calculates Eigen Centrality Z-Score, Shortest Path Z-score and
Betweenness Z-Score of input positions. Find the mutation type
(1-Deadly, 0-Not Deadly), find the Gene Name from the PDB file and if
there are multiple genes annotated for the input PDB file, ask user to
choose the input gene, gets GnomAD Synonymous Z-Score, Non-Synonymous
Z-Score and PLoF Score, Genic Intolerance Score, gets BLOSUM62 score of
the amino acid change, finds the KEGG Pathway number which contains the
input gene. Make prediction based on an Adaboost model and classifies
the mutation as *Disease Causing* or *Neutral*.
