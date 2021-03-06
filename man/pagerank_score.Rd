% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pagerank_score.R
\name{pagerank_score}
\alias{pagerank_score}
\title{PageRank Score}
\usage{
pagerank_score(edge_list, filtered_info_df)
}
\arguments{
\item{edge_list}{list contains separate edge data.frames for each chain}

\item{filtered_info_df}{input data.frame which contain only one PDB entries}
}
\value{
PageRank Z-Scores of input position
}
\description{
\code{pagerank_score} function assigns an importance score to every node via
\code{page_rank} function of igraph package. Calculates and returns the Z-Scores of the input positions.
}
