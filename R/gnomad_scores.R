#' GnomAD Scores
#'
#' @param PDB_ID PDB ID
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return list that contains a data.frame containing GnomAD scores of input mutations and
#' the gene name
#'
#' @export
#'

gnomad_scores <- function(PDB_ID, filtered_info_df){

  if(length(unique(filtered_info_df$Gene_Name)) == 1 & "" %in% unique(filtered_info_df$Gene_Name)){

  res <- biomart_data[which(biomart_data$pdb == PDB_ID),]

  if(length(res$external_gene_name) > 1){
    gene_name <- readline(prompt = paste0("There are multiple gene names for PDB ", PDB_ID,
                                          ", please write the input gene name (",
                                          paste0(res$external_gene_name, collapse = ", "),
                                          "): "))

    if(!(gene_name %in% res$external_gene_name)){

      stop(paste0("Input gene name ", gene_name,
                  " not found in query genes (",
                  paste0(res$external_gene_name, collapse = ", "), ")"))
    }

    idx <- which(gnomad_data$gene == gene_name)

    if(length(idx) == 0){
      message(crayon::white(paste0("\n", "gnomAD scores of the gene ", gene_name,
                                   " couldn't find, will be removed from the query", "\n")))

      filtered_info_df$syn_z <- NA

      filtered_info_df$mis_z <- NA

      filtered_info_df$pLI <- NA

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(filtered_info_df, gene_name))
    }
    else{

      filtered_info_df$syn_z <- max(as.numeric(gnomad_data$syn_z[idx]))

      filtered_info_df$mis_z <- max(as.numeric(gnomad_data$mis_z[idx]))

      filtered_info_df$pLI <- max(as.numeric(gnomad_data$pLI[idx]))

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(filtered_info_df, gene_name))
    }

  }
  else if(length(res$external_gene_name) == 0){
    gene_name <- readline(prompt = paste0("No gene name found for the PDB ", PDB_ID,
                                          ", please write the input gene name: "))

    idx <- which(gnomad_data$gene == gene_name)

    if(length(idx) == 0){
      message(crayon::white(paste0("\n", "gnomAD scores of the gene ", gene_name,
                                   " couldn't find, will be removed from the query", "\n")))

      filtered_info_df$syn_z <- NA

      filtered_info_df$mis_z <- NA

      filtered_info_df$pLI <- NA

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(filtered_info_df, gene_name))
    }
    else{

      filtered_info_df$syn_z <- max(as.numeric(gnomad_data$syn_z[idx]))

      filtered_info_df$mis_z <- max(as.numeric(gnomad_data$mis_z[idx]))

      filtered_info_df$pLI <- max(as.numeric(gnomad_data$pLI[idx]))

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(filtered_info_df, gene_name))
    }
  }
  else{
    idx <- which(gnomad_data$gene == res$external_gene_name)

    filtered_info_df$syn_z <- max(as.numeric(gnomad_data$syn_z[idx]))

    filtered_info_df$mis_z <- max(as.numeric(gnomad_data$mis_z[idx]))

    filtered_info_df$pLI <- max(as.numeric(gnomad_data$pLI[idx]))

    message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

    return(list(filtered_info_df, gnomad_data$gene[idx]))
  }
  }
  else{
    if(length(unique(filtered_info_df$Gene_Name)) != 1){
      message(crayon::white(paste0("\n", "Gene name input should be same for the same PDB ID (", unique(filtered_info_df$PDB_ID), "), inputs will be removed from the query", "\n")))

      filtered_info_df$syn_z <- NA

      filtered_info_df$mis_z <- NA

      filtered_info_df$pLI <- NA

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(filtered_info_df, "no_name"))
    }
    else{
      idx <- which(gnomad_data$gene == unique(filtered_info_df$Gene_Name))

      if(length(idx) == 0){
        message(crayon::white(paste0("\n", "gnomAD scores of the gene ", unique(filtered_info_df$Gene_Name),
                                     " couldn't find, will be removed from the query", "\n")))

        filtered_info_df$syn_z <- NA

        filtered_info_df$mis_z <- NA

        filtered_info_df$pLI <- NA

        message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

        return(list(filtered_info_df, unique(filtered_info_df$Gene_Name)))
      }
      else{

        filtered_info_df$syn_z <- max(as.numeric(gnomad_data$syn_z[idx]))

        filtered_info_df$mis_z <- max(as.numeric(gnomad_data$mis_z[idx]))

        filtered_info_df$pLI <- max(as.numeric(gnomad_data$pLI[idx]))

        message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

        return(list(filtered_info_df, unique(filtered_info_df$Gene_Name)))
      }
    }
  }
}


