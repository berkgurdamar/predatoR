#' PDB Reading/Downloading
#'
#' Download and read input PDB via \code{get.pdb} and \code{read.pdb} function of bio3d package
#' and filter the PDB file for only containing atoms.
#'
#' @param PDB_ID PDB ID
#' @param internal_PDB_path PDB file path (default = FALSE)
#'
#' @return Matrix that contains all the atoms in the PDB structure
#' @export
#'

read_PDB <- function(PDB_ID, internal_PDB_path = NULL){

  if(is.null(internal_PDB_path) == TRUE){

    try <- tryCatch({
      pdb <- suppressMessages(bio3d::read.pdb(bio3d::get.pdb(PDB_ID, URL=TRUE)))

      atom_matrix <- pdb$atom
      atom_matrix <- atom_matrix[which(atom_matrix$type == "ATOM"),]
      atom_matrix <- atom_matrix[which(atom_matrix$resid != "HOH"),]

      message(crayon::white(paste0("STEP:","\n" ,"Reading PDB:", "\t\t\t", "DONE")))

      return(atom_matrix)

    },
    error=function(x){
      message(paste0("\n", "Couldn't download the PDB file ", PDB_ID, ", it will be removed from the query", "\n"))
      return(NA)
    },
    warning=function(x) {
      message(paste0("\n", "Couldn't download the PDB file ", PDB_ID, ", it will be removed from the query", "\n"))
      return(NA)
    }
    )
  }else{
    pdb <- suppressMessages(bio3d::read.pdb(paste0(internal_PDB_path, PDB_ID, ".pdb")))

    atom_matrix <- pdb$atom
    atom_matrix <- atom_matrix[which(atom_matrix$type == "ATOM"),]
    atom_matrix <- atom_matrix[which(atom_matrix$resid != "HOH"),]

    message(crayon::white(paste0("STEP:","\n" ,"Reading PDB:", "\t\t\t", "DONE")))

    return(atom_matrix)

  }
}
