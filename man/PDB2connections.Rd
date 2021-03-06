% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDB2connections.R
\name{PDB2connections}
\alias{PDB2connections}
\title{Creating edge list from atom matrix}
\usage{
PDB2connections(
  atom_matrix,
  filtered_info_df,
  n_threads = NULL,
  single_run = TRUE
)
}
\arguments{
\item{atom_matrix}{matrix that contains all the atoms in the PDB file}

\item{filtered_info_df}{input data.frame which contain only one PDB entries}

\item{n_threads}{number of threads (default = NULL)}

\item{single_run}{should be set as TRUE when using \code{PDB2connections} function alone (default = TRUE)}
}
\value{
list contains separate edge data.frames for each chain
}
\description{
Creates all the edges from the input atom matrix created by \code{PDB_read} .
}
\details{
This function calculates the distances between all the atoms in the atom matrix and create edges
if the distance between two atoms is less than or equal to 7 Angstrom.
}
