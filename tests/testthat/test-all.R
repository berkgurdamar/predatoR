library(testthat)
library(predatoR)


# predatoR ----------------------------------------------------------------



info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))


test_that("Check input type", {
  expect_error(predatoR(info_df = as.matrix(info_df), gene_name_info = T),
               "Input should be a data.frame.")
})

test_that("Check input column number", {
  expect_error(predatoR(info_df = info_df, gene_name_info = F),
                 "Input data.frame should contain 5 columns; PDB_ID, Chain, Position, Orig_AA and Mut_AA respectively.")
})

test_that("Check input column number", {
  expect_error(predatoR(info_df = info_df[,1:5], gene_name_info = T),
               "Input data.frame should contain 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name respectively.")
})

test_that("Check output class and thread input", {

  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "HBB"),
                                 c("4ONL", "A", 36, "GLU", "GLY", "UBE2V2"),
                                 c("4ONL", "A", 78, "PRO", "GLN", "UBE2V2")))

  expect_true(is.data.frame(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2)))
})


test_that("Check false pdb input", {

  info_df <- as.data.frame(rbind(c("2DNNN2", "B", 1, "VAL", "ALA", "HBB")))

  expect_error(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2),
               "There is no input for prediction")
})


test_that("Check no gene name input", {

  info_df <- as.data.frame(rbind(c("2DNNN2", "B", 1, "VAL", "ALA")))

  expect_error(predatoR(info_df = info_df, gene_name_info = F, n_threads = 2),
               "There is no input for prediction")
})


### fix regex
info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "asdad", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

test_that("Check input amino acid names", {
  expect_error(predatoR(info_df = info_df, gene_name_info = T))
})

##

info_df <- as.data.frame(rbind(c("2DN2", "B", 10000, "VAL", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

test_that("Check if residue included", {
  expect_message(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2),
                 "Residue [1-9]\\d* is not included in the PDB structure, it will be removed from the query")
})


info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "GLU", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

test_that("Check if residue-amino acid matching", {
  expect_message(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2),
                 "Residue [1-9]\\d* is not [a-zA-Z]+ in the PDB structure, it will be removed from the query"
  )
})


info_df <- as.data.frame(rbind(c("2DN2", "B", 5, "VAL", "ALA", "HBB"),
                               c("2DN2", "B", 10, "GLU", "ALA", "HBB")))


test_that("Check input amino acid names", {
  expect_error(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2),
               "There is no input for prediction")
})



test_that("Check multiple amino acid for one residue message", {

  info_df <- as.data.frame(rbind(c("2DN2_N", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DN2_N", "B", 6, "GLU", "ALA", "HBB")))
  PDB_ID <- "2DN2"
  tmp_dir <- tempdir()
  tmp_path <- file.path(tmp_dir, paste0(PDB_ID, ".pdb"))
  download.file(bio3d::get.pdb(PDB_ID, URL=TRUE), tmp_path)
  pdb_file <- bio3d::read.pdb(tmp_path)

  pdb_file$atom$resid[1070] <- "ASP"

  bio3d::write.pdb(pdb_file, file.path(tmp_dir, paste0(PDB_ID, "_N.pdb")))

  expect_message(predatoR(info_df = info_df, gene_name_info = T, PDB_path = tmp_dir, n_threads = 2),
                 "There are multiple amino acids for residue [1-9]\\d* in the PDB file, it will be removed from the query")

})


# read_PDB ----------------------------------------------------------------

test_that("Check local PDB read", {
  PDB_ID <- "2DN2"
  tmp_dir <- tempdir()
  tmp_path <- file.path(tmp_dir, paste0(PDB_ID, ".pdb"))
  download.file(bio3d::get.pdb(PDB_ID, URL=TRUE), tmp_path)

  expect_true(is.data.frame(read_PDB(PDB_ID = PDB_ID, PDB_path = tmp_dir)))

})


test_that("Check wrong PDB name", {
  expect_message(read_PDB("ACCSCS"),
                 "Couldn't download the PDB file [a-zA-Z]+, it will be removed from the query"
  )
})


test_that("Check wrong path", {
  expect_message(read_PDB("2DN2", PDB_path = "path/"),
                 "Couldn't find the PDB file [A-Za-z0-9_-]+ in the given path, it will be removed from the query"
  )
})

test_that("Check output class", {
  expect_true(is.data.frame(read_PDB("2DN2")))
})



# PDB2connections ---------------------------------------------------------

atom_matrix <- read_PDB("2DN2")

info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

connections_df <- PDB2connections(atom_matrix, info_df, single_run = TRUE, n_threads = 2)

test_that("Check output class an thread input", {

  expect_true(is.list(connections_df))
})


test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(PDB2connections(atom_matrix, info_df, single_run = TRUE),
               "filtered_info_df should contain only one PDB entries")
  expect_error(PDB2connections(atom_matrix, info_df, single_run = FALSE),
               "filtered_info_df should contain only one PDB entries")
})




# eigen_centrality_score --------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(eigen_centrality_score(connections_df, info_df)))
  expect_message(eigen_centrality_score(connections_df, info_df),
                 "Eigen Centrality Score:		DONE")
})

test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(eigen_centrality_score(connections_df, info_df),
               "filtered_info_df should contain only one PDB entries")
})




# shorteset_path_score ----------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(shorteset_path_score(connections_df, info_df)))
  expect_message(shorteset_path_score(connections_df, info_df),
                 "Shortest Path Score:		DONE")
})

test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(shorteset_path_score(connections_df, info_df),
               "filtered_info_df should contain only one PDB entries")
})


# betweenness_score -------------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(betweenness_score(connections_df, info_df)))
  expect_message(betweenness_score(connections_df, info_df),
                 "Betweenness Score:		DONE")
})

test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(betweenness_score(connections_df, info_df),
               "filtered_info_df should contain only one PDB entries")
})



# gnomad_scores -----------------------------------------------------------


test_that("Check output class", {

  expect_true(is.data.frame(gnomad_scores(info_df)))
  expect_message(gnomad_scores(info_df),
                 "GNOMAD Scores:			DONE")
})

test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DNN", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(gnomad_scores(info_df),
               "filtered_info_df should contain only one PDB entries")
})


test_that("Check column number error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA"),
                                 c("2DN2", "B", 6, "GLU", "ALA")))

  expect_error(gnomad_scores(info_df),
               "Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
})



test_that("Check wrong gene name message", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "BRCA"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "BRCA")))

  expect_message(gnomad_scores(info_df))

})


test_that("Check multiple gene info", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "AQP1"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "AQP1")))

  expect_true(is.data.frame(gnomad_scores(info_df)))

})


test_that("Check no gene input returning multiple gene", {
  info_df <- as.data.frame(rbind(c("1A01", "A", 1, "VAL", "ALA", ""),
                                 c("1A01", "A", 6, "GLU", "ALA", "")))

  expect_true(is.data.frame(gnomad_scores(info_df)))
})

test_that("Check no gene input", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", ""),
                                 c("2DN2", "B", 6, "GLU", "ALA", "")))

  expect_true(is.data.frame(gnomad_scores(info_df)))
})

test_that("Check wrong pdb input returning no gene", {
  info_df <- as.data.frame(rbind(c("ASDA", "B", 1, "VAL", "ALA", ""),
                                 c("ASDA", "B", 6, "GLU", "ALA", "")))

  expect_true(is.data.frame(gnomad_scores(info_df)))
})


test_that("Check no gene info in gnomad", {
  info_df <- as.data.frame(rbind(c("1ND5", "B", 1, "VAL", "ALA", ""),
                                 c("1ND5", "B", 6, "GLU", "ALA", "")))

  expect_true(is.data.frame(gnomad_scores(info_df)))
})


# BLOSUM62_score ----------------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(BLOSUM62_score(info_df)))
  expect_message(BLOSUM62_score(info_df),
                 "BLOSUM62 Score:			DONE")
})


# KEGG_pathway_number -----------------------------------------------------

test_that("Check output class", {

  expect_true(is.data.frame(KEGG_pathway_number(info_df)))
  expect_message(KEGG_pathway_number(info_df),
                 "KEGG Pathway Number:		DONE")
})



test_that("Check column number error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA"),
                                 c("2DN2", "B", 6, "GLU", "ALA")))

  expect_error(KEGG_pathway_number(info_df),
               "Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
})


# genic_intolerance -------------------------------------------------------

test_that("Check output class", {

  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

  expect_true(is.data.frame(genic_intolerance(info_df)))
  expect_message(genic_intolerance(info_df),
                 "Genic Intolerance Score:	DONE")
})



test_that("Check input error message", {

  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "BRCA"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "BRCA")))

  expect_message(genic_intolerance(info_df))
})



test_that("Check column number error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA"),
                                 c("2DN2", "B", 6, "GLU", "ALA")))

  expect_error(genic_intolerance(info_df),
               "Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
})

# clique_score ------------------------------------------------------------


atom_matrix <- read_PDB("2DN2")

info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                               c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

connections_df <- PDB2connections(atom_matrix, info_df, single_run = TRUE, n_threads = 2)



test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(clique_score(connections_df, info_df, single_run = TRUE),
               "filtered_info_df should contain only one PDB entries")
  expect_error(clique_score(connections_df, info_df, single_run = FALSE),
               "filtered_info_df should contain only one PDB entries")
})


test_that("Check output class", {

  expect_true(is.numeric(clique_score(connections_df, info_df, single_run = TRUE, n_threads = 2)))
})




# degree_score ------------------------------------------------------------




test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(degree_score(connections_df, info_df),
               "filtered_info_df should contain only one PDB entries")
})


test_that("Check output class", {

  expect_true(is.numeric(degree_score(connections_df, info_df)))
})





# pagerank_score ----------------------------------------------------------


test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(pagerank_score(connections_df, info_df),
               "filtered_info_df should contain only one PDB entries")
})


test_that("Check output class", {

  expect_true(is.numeric(pagerank_score(connections_df, info_df)))
})



# impact_prediction -------------------------------------------------------

final_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB", "-0.3319392", "-0.9048452", "0.1615106", "-0.04987094",
                                  "-0.8371797", "-0.1839022", "-3.7953", "-0.2113", "1.2274e-09", "0", "5", "-0.01"),
                                c("2DN2", "B", 6, "GLU", "ALA", "HBB", "-0.6075531", "-1.3857613", "1.8748888", "-0.86584380",
                                  "-0.3499602", "-0.3829808", "-3.7953", "-0.2113", "1.2274e-09", "-1", "5", "-0.01")))

colnames(final_df) <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name", "degree_z_score",
                        "eigen_z_score", "shortest_path_z", "betweenness_scores_z", "clique_z_score", "pagerank_z_score",
                        "syn_z", "mis_z", "pLI", "blosum62_scores", "kegg_pathway_number", "genic_intolerance")


test_that("Check output class", {

  expect_true(is.data.frame(impact_prediction(final_df)))

})
