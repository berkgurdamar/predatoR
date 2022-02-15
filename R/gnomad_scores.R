#' GnomAD Scores
#'
#' Function gets input gene's gnomAD Synonymous Z-Score, Non-Synonymous Z-Score, and PLoF Score
#' from data retrieved from [gnomAD](https://gnomad.broadinstitute.org/). If there is no gene name information in the input,
#' \code{gnomad_scores} finds the genes that related with input PDB ID by using a dataset retrieved from [Ensembl](https://www.ensembl.org/).
#' If there are multiple genes annotated for same PDB, \code{gnomad_scores} function gets the gene that has maximum Synonymous Z-Score, Non-Synonymous Z-Score, and PLoF Score.
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return gnomAD scores added input data.frame
#'
#' @export
#'

gnomad_scores <- function(filtered_info_df){

  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }

  if(length(unique(filtered_info_df$Gene_Name)) == 2 & "" %in% unique(filtered_info_df$Gene_Name)){
    gene_names <- unique(filtered_info_df$Gene_Name[filtered_info_df$Gene_Name != ""])
    filtered_info_df$Gene_Name <- gene_names

  }


  no_name <- c()
  no_gnomad <- c()
  for(i in 1:length(filtered_info_df$Gene_Name)){

    if(filtered_info_df$Gene_Name[i] == ""){

      filtered_info_df$pdb_chain_info <- paste0(filtered_info_df$PDB_ID, ".", filtered_info_df$Chain)

      chain_info <- unique(filtered_info_df$pdb_chain_info)

      res <- pdb2gene[which(pdb2gene$PDB.Chain %in% chain_info),]

      if(length(unique(res$Gene)) == 0){

        no_name <- c(no_name, paste0(unique(filtered_info_df$pdb_chain_info), collapse = ", "))

        filtered_info_df$Gene_Name[i] <- "no_name"

        filtered_info_df$syn_z[i] <- NA

        filtered_info_df$mis_z[i] <- NA

        filtered_info_df$pLI[i] <- NA

      }
      else if(length(unique(res$Gene)) == 1){

        filtered_gnomad_data <- gnomad_data[gnomad_data$gene == unique(res$Gene),]

        if(nrow(filtered_gnomad_data) == 0){

          no_gnomad <- c(no_gnomad, unique(res$Gene))

          filtered_info_df$Gene_Name[i] <- unique(res$Gene)

          filtered_info_df$syn_z[i] <- NA

          filtered_info_df$mis_z[i] <- NA

          filtered_info_df$pLI[i] <- NA
        }
        else{

          idx <- table(sapply(2:4, function(x) which.max(filtered_gnomad_data[,x])))
          filtered_gnomad_data <- filtered_gnomad_data[as.numeric(names(idx)[which.max(idx)]),]

          filtered_info_df$Gene_Name[i] <- filtered_gnomad_data$gene

          filtered_info_df$syn_z[i] <- filtered_gnomad_data$syn_z

          filtered_info_df$mis_z[i] <- filtered_gnomad_data$mis_z

          filtered_info_df$pLI[i] <- filtered_gnomad_data$pLI
        }
      }
      else{

        filtered_gnomad_data <- gnomad_data[gnomad_data$gene %in% unique(res$Gene),]

        if(nrow(filtered_gnomad_data) == 0){

          no_gnomad <- c(no_gnomad, paste(unique(res$Gene), collapse = ", "))

          filtered_info_df$Gene_Name[i] <- paste(unique(res$Gene), collapse = ", ")

          filtered_info_df$syn_z[i] <- NA

          filtered_info_df$mis_z[i] <- NA

          filtered_info_df$pLI[i] <- NA
        }
        else{

        idx <- table(sapply(2:4, function(x) which.max(filtered_gnomad_data[,x])))
        filtered_gnomad_data <- filtered_gnomad_data[as.numeric(names(idx)[which.max(idx)]),]

        filtered_info_df$Gene_Name[i] <- filtered_gnomad_data$gene

        filtered_info_df$syn_z[i] <- filtered_gnomad_data$syn_z

        filtered_info_df$mis_z[i] <- filtered_gnomad_data$mis_z

        filtered_info_df$pLI[i] <- filtered_gnomad_data$pLI
        }
      }

      filtered_info_df <- filtered_info_df[-(which(colnames(filtered_info_df) == "pdb_chain_info"))]


    }

    # gene name included
    else{

      filtered_gnomad_data <- gnomad_data[gnomad_data$gene == filtered_info_df$Gene_Name[i],]

      if(nrow(filtered_gnomad_data) == 0){

        no_gnomad <- c(no_gnomad, filtered_info_df$Gene_Name[i])

        filtered_info_df$syn_z[i] <- NA

        filtered_info_df$mis_z[i] <- NA

        filtered_info_df$pLI[i] <- NA
      }

      else if(nrow(filtered_gnomad_data) == 1){

        filtered_info_df$syn_z[i] <- filtered_gnomad_data$syn_z

        filtered_info_df$mis_z[i] <- filtered_gnomad_data$mis_z

        filtered_info_df$pLI[i] <- filtered_gnomad_data$pLI

      }
      else{

        idx <- table(sapply(2:4, function(x) which.max(filtered_gnomad_data[,x])))
        filtered_gnomad_data <- filtered_gnomad_data[as.numeric(names(idx)[which.max(idx)]),]

        filtered_info_df$syn_z[i] <- filtered_gnomad_data$syn_z

        filtered_info_df$mis_z[i] <- filtered_gnomad_data$mis_z

        filtered_info_df$pLI[i] <- filtered_gnomad_data$pLI
      }
    }
  }

  if(length(unique(no_name)) > 0){
    message(crayon::white(paste0("\n", "Gene name(s) couldn't find for ", paste(unique(no_name), collapse = ", "), " (PDB.Chain), will be removed from the query", "\n")))
  }

  if(length(unique(no_gnomad)) > 0){
    message(crayon::white(paste0("\n", "gnomAD scores of the gene(s) ", paste(unique(no_gnomad), collapse = ", "), " couldn't find, will be removed from the query", "\n")))
  }

  message(crayon::white(paste0("GNOMAD Scores:", "\t\t\t", "DONE")))

  return(filtered_info_df)
}


