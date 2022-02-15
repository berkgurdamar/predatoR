#' PSSM Scores
#'
#'
#' @param filtered_info_df input data.frame which contain PDB entries
#' @param pssm_path folder path for saving PSSM matrices for fastening the further analysis
#'
#' @return PSSM scores of input positions
#' @export
#'


PSSM_Scores <- function(filtered_info_df, pssm_path = NULL){

  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  if(is.null(pssm_path)){
    pssm_path <- file.path(.libPaths(), "predatoR", "pssm_matrix.RDS")
    if(file.exists(pssm_path)){
      pssm_matrix <- readRDS(pssm_path)
    }
    else{
      pssm_matrix <- list()
    }
  }else{
    pssm_path <- file.path(pssm_path, "pssm_matrix.RDS")
    if(file.exists(pssm_path)){
      pssm_matrix <- readRDS(pssm_path)
    }
    else{
      pssm_matrix <- list()
    }
  }

  get_seq <- purrr::possibly(bio3d::get.seq, otherwise = "no_sequence")
  hmmer_new <- purrr::possibly(bio3d::hmmer, otherwise = "no_hmmer")

  dbpath <- tempfile("tempdb", fileext = ".fasta")

  filtered_info_df$Position <- as.numeric(filtered_info_df$Position)
  filtered_info_df$pssm <- NA

  for(i in 1:nrow(filtered_info_df)){

    print(i)

    if(paste0(filtered_info_df$PDB_ID[i], "_",
              filtered_info_df$Chain[i]) %in% names(pssm_matrix)){

      pssm <- pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
                                  filtered_info_df$Chain[i])]]

      if(!is.data.frame(pssm)){
        filtered_info_df$pssm[i] <- NA
        next
      }

      score <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
        pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]

      if(length(score) == 0){

        filtered_info_df$pssm[i] <- NA

      }
      else{

        filtered_info_df$pssm[i] <- score

      }


    }else{

      invisible(capture.output(seq <- suppressWarnings(get_seq(paste0(filtered_info_df$PDB_ID[i], "_", filtered_info_df$Chain[i])))))

      if("no_sequence" %in% seq){
        pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
                            filtered_info_df$Chain[i])]] <- NA

        filtered_info_df$pssm[i] <- NA
        next
      }
      hmm <- hmmer_new(seq, type="hmmscan", db="pfam")

      if("no_hmmer" %in% hmm){
        pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
                            filtered_info_df$Chain[i])]] <- NA

        filtered_info_df$pssm[i] <- NA
        next
      }

      aln <- bio3d::pfam(hmm$hit.tbl$acc[which.min(hmm$hit.tbl$evalue)])

      sequences <- sapply(1:nrow(aln$ali), function(x) gsub("-", "", paste0(aln$ali[x,], collapse = "")))

      writeLines(paste(paste0(">", aln$id), "\n", sequences, collapse = "\n"), dbpath)

      pssm <- as.data.frame(protr::extractPSSM(paste0(seq$ali, collapse = ""), database.path = dbpath))

      # pdb_sequences[[paste0(filtered_info_df$PDB_ID[i], "_", filtered_info_df$Chain[i])]] <- paste0(seq$ali, collapse = "")

      score <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
        pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]

      if(length(score) == 0){

        filtered_info_df$pssm[i] <- NA

      }
      else{

        filtered_info_df$pssm[i] <- score

      }

      pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
                          filtered_info_df$Chain[i])]] <- pssm

    }

  }
    saveRDS(pssm_matrix, pssm_path)

    message(crayon::white(paste0("PSSM Scores:", "\t\t\t", "DONE")))

    return(filtered_info_df)


}

# PSSM_Scores <- function(filtered_info_df, pssm_path = NULL){
#
#   if(ncol(filtered_info_df) < 6){
#     stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
#   }
#   else{
#     colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
#   }
#
#   if(is.null(pssm_path)){
#     pssm_path <- file.path(.libPaths(), "predatoR", "pssm_matrix.RDS")
#     if(file.exists(pssm_path)){
#       pssm_matrix <- readRDS(pssm_path)
#     }
#     else{
#       pssm_matrix <- list()
#     }
#   }else{
#     pssm_path <- file.path(pssm_path, "pssm_matrix.RDS")
#     if(file.exists(pssm_path)){
#       pssm_matrix <- readRDS(pssm_path)
#     }
#     else{
#       pssm_matrix <- list()
#     }
#   }
#
#   get_seq <- purrr::possibly(bio3d::get.seq, otherwise = "no_sequence")
#
#   dbpath <- tempfile("tempdb", fileext = ".fasta")
#
#   filtered_info_df$Position <- as.numeric(filtered_info_df$Position)
#
#   for(i in 1:nrow(filtered_info_df)){
#
#     print(i)
#     tmp <- which(pfam_db$PDB_Chain %in% paste0(filtered_info_df$PDB_ID[i], "_", filtered_info_df$Chain[i]))
#
#     if(length(tmp) == 0){
#
#       filtered_info_df$pssm[i] <- NA
#
#     }else{
#       filtered_pfam <- pfam_db[tmp,]
#       idx <- filtered_pfam$PFAM_ACCESSION[filtered_pfam$PDB_START <= filtered_info_df$Position[i] & filtered_pfam$PDB_END >= filtered_info_df$Position[i]]
#
#       if(length(idx) == 1){
#
#         if(paste0(filtered_info_df$PDB_ID[i], "_",
#                   filtered_info_df$Chain[i], "_",
#                   idx) %in% names(pssm_matrix)){
#
#           pssm <- pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
#                                       filtered_info_df$Chain[i], "_",
#                                       idx)]]
#
#           score <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
#             pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]
#
#           if(length(score) == 0){
#
#             filtered_info_df$pssm[i] <- NA
#
#           }
#           else{
#
#             filtered_info_df$pssm[i] <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
#               pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]
#
#           }
#
#
#         }else{
#
#           invisible(capture.output(seq <- suppressWarnings(get_seq(paste0(filtered_info_df$PDB_ID[i], "_", filtered_info_df$Chain[i])))))
#
#           if("no_sequence" %in% seq){
#             next
#           }
#
#           aln <- bio3d::pfam(idx)
#
#           sequences <- sapply(1:nrow(aln$ali), function(x) gsub("-", "", paste0(aln$ali[x,], collapse = "")))
#
#           writeLines(paste(paste0(">", aln$id), "\n", sequences, collapse = "\n"), dbpath)
#
#           pssm <- as.data.frame(protr::extractPSSM(paste0(seq$ali, collapse = ""), database.path = dbpath))
#
#           score <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
#             pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]
#
#           if(length(score) == 0){
#
#             filtered_info_df$pssm[i] <- NA
#
#           }
#           else{
#
#
#             filtered_info_df$pssm[i] <- pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Mut_AA[i])),paste0("V", filtered_info_df$Position[i])] -
#               pssm[seqinr::a(stringr::str_to_title(filtered_info_df$Orig_AA[i])),paste0("V", filtered_info_df$Position[i])]
#
#           }
#
#           pssm_matrix[[paste0(filtered_info_df$PDB_ID[i], "_",
#                               filtered_info_df$Chain[i], "_",
#                               idx)]] <- pssm
#
#         }
#
#       }else{
#
#         filtered_info_df$pssm[i] <- NA
#       }
#
#     }
#
#   }
#   saveRDS(pssm_matrix, pssm_path)
#
#   message(crayon::white(paste0("PSSM Scores:", "\t\t\t", "DONE")))
#
#   return(filtered_info_df)
#
#
# }

