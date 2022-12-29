#' PDB Reading/Downloading
#'
#' Download/Read input PDB via \code{get.pdb} and \code{read.pdb} function of bio3d package
#' and filter the PDB file for only containing atoms.
#'
#' @param PDB_ID PDB ID
#' @param PDB_path PDB file path (default = NULL)
#' @param network_approach network building approach; "ca" (default) for using ca atoms only or "all" for using all atoms
#'
#' @return Matrix that contains all the atoms in the PDB structure
#'
#' @export
#'

read_PDB <- function(PDB_ID, PDB_path = NULL, network_approach = "ca"){

  if(!network_approach %in% c("all", "ca")){
    stop("Network approach needs to be 'all' or 'ca'")

  }

  if(is.null(PDB_path)){

    try <- tryCatch({
      pdb <- suppressMessages(bio3d::read.pdb(bio3d::get.pdb(PDB_ID, URL=TRUE)))

      atom_matrix <- pdb$atom
      atom_matrix <- atom_matrix[which(atom_matrix$type == "ATOM"),]
      atom_matrix <- atom_matrix[which(atom_matrix$resid != "HOH"),]

      if(network_approach == "ca"){
        atom_matrix <- atom_matrix[which(atom_matrix$elety == "CA"),]
      }

      message(crayon::white(paste0("STEP:","\n" ,"Reading PDB:", "\t\t\t", "DONE")))

      return(atom_matrix)

    },

    warning=function(x){

      message(crayon::white(paste0("\n", "Couldn't download the PDB file ",
                                   PDB_ID,
                                   ", it will be removed from the query", "\n")))
      return(NA)

    }
    )
  } else {

    f_path <- file.path(PDB_path, paste0(PDB_ID, ".pdb"))

    if (file.exists(f_path)) {

      pdb <- suppressMessages(bio3d::read.pdb(f_path))

      atom_matrix <- pdb$atom
      atom_matrix <- atom_matrix[which(atom_matrix$type == "ATOM"),]
      atom_matrix <- atom_matrix[which(atom_matrix$resid != "HOH"),]

      if(network_approach == "ca"){
        atom_matrix <- atom_matrix[which(atom_matrix$elety == "CA"),]
      }

      message(crayon::white(paste0("STEP:","\n" ,"Reading PDB:", "\t\t\t", "DONE")))

      return(as.data.frame(atom_matrix))

    } else {
      message(crayon::white(paste0("\n", "Couldn't find the PDB file ",
                                   PDB_ID,
                                   " in the given path, it will be removed from the query", "\n")))
      return(NA)

    }
  }
}
