# Sys.setenv(R_TESTS="")

library(testthat)
library(predatoR)

# chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#
# if (nzchar(chk) && chk == "TRUE") {
#   # use 2 cores in CRAN/Travis/AppVeyor
#   num_workers <- 2L
# } else {
#   # use all cores in devtools::test()
#   num_workers <- parallel::detectCores() - 1
# }
# predatoR ----------------------------------------------------------------

test_that("Check wrong gene name message", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "BRCA"),
                                 c("4ONL", "A", 36, "GLU", "GLY", "UBE2V2")))

  expect_message(predatoR(info_df = info_df, gene_name_info = T, n_threads = 2),
                 "Gene name input should be same for the same PDB ID [()A-Za-z0-9_-]+, inputs will be removed from the query")

})

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
                                 c("2DN2", "B", 6, "GLU", "ALA", "HBB")))

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
  expect_message(PDB2connections(atom_matrix, info_df, single_run = TRUE, n_threads = 2),
                 "List of Edges:			DONE")
})


test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DDD", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(PDB2connections(atom_matrix, info_df, single_run = TRUE),
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
                 "Number of Shortest Paths:	DONE")
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

  expect_true(is.list(gnomad_scores("2DN2", info_df)))
  expect_message(gnomad_scores("2DN2", info_df),
                 "GNOMAD Information:		DONE")
})

test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DNN", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(gnomad_scores("2DN2", info_df),
               "filtered_info_df should contain only one PDB entries")
})


test_that("Check column name error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA"),
                                 c("2DN2", "B", 6, "GLU", "ALA")))

  expect_error(gnomad_scores("2DN2", info_df),
               "Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
})



test_that("Check multiple PDB error", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DNN", "B", 6, "GLU", "ALA", "HBB")))

  expect_error(gnomad_scores("2DN2", info_df),
               "filtered_info_df should contain only one PDB entries")
})

test_that("Check wrong gene name message", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "ASDSADASD"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "ASDSADASD")))

  expect_message(gnomad_scores("2DN2", info_df),
               "gnomAD scores of the gene [A-Za-z0-9_-]+ couldn't find, will be removed from the query")

})



test_that("Check multiple gene name message", {
  info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
                                 c("2DN2", "B", 6, "GLU", "ALA", "BRCA")))

  expect_message(gnomad_scores("2DN2", info_df),
                 "Gene name input should be same for the same PDB ID [()A-Za-z0-9_-]+, inputs will be removed from the query")

})


# find pdb for this
# test_that("Check multiple amino acid for a residue message", {
#   info_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB"),
#                                  c("2DN2", "B", 6, "GLU", "ALA", "BRCA")))
#
#   expect_message(gnomad_scores("2DN2", info_df),
#                  "Gene name input should be same for the same PDB ID [()A-Za-z0-9_-]+, inputs will be removed from the query")
#   expect_true(gnomad_scores("2DN2", info_df)[[2]] == "no_name")
#
# })



# BLOSUM62_score ----------------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(BLOSUM62_score(info_df)))
  expect_message(BLOSUM62_score(info_df),
                 "BLOSUM62 Score:			DONE")
})


# KEGG_pathway_number -----------------------------------------------------

test_that("Check output class", {

  expect_true(is.numeric(KEGG_pathway_number("HBB")))
  expect_message(KEGG_pathway_number("HBB"),
                 "Number of KEGG Pathways:	DONE")
})




# genic_intolerance -------------------------------------------------------

test_that("Check output class", {

  expect_true(is.character(genic_intolerance("HBB")))
  expect_message(genic_intolerance("HBB"),
                 "Genic Intolerance Score:	DONE")
})



test_that("Check input error message", {

  expect_message(genic_intolerance("ASDASDASD"),
                 "Genic Intolerance score of the gene [A-Za-z0-9_-]+ couldn't find, will be removed from the query")
})


test_that("Check input error output type", {

  expect_true(is.na(genic_intolerance("ASDASDASD")))
})




# impact_prediction -------------------------------------------------------

final_df <- as.data.frame(rbind(c("2DN2", "B", 1, "VAL", "ALA", "HBB", "-0.7491699", "0.1615106", "-0.04987094",
                                  "-3.7953", "-0.2113", "1.2274e-09", "0", "5", "-0.01"),
                                c("2DN2", "B", 6, "GLU", "ALA", "HBB", "-1.2004364", "1.8748888", "-0.86584380",
                                  "-3.7953", "-0.2113", "1.2274e-09", "-1", "5", "-0.01")))

colnames(final_df) <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name",
                        "eigen_z_score", "shortest_path_z", "betwenness_scores_z", "syn_z", "mis_z", "pLI",
                        "blosum62_scores", "kegg_pathway_number", "genic_intolerance")


test_that("Check output class", {

  expect_true(is.data.frame(impact_prediction(final_df)))

})
