#' PDB Downloading and Reading
#'
#' Download and read input PDB via \code{get.pdb} and \code{read.pdb} function of bio3d package
#' and filter the PDB file for only containing atoms.
#'
#' @param PDB_ID PDB ID
#'
#' @return Matrix that contains all the atoms in the PDB structure
#' @export
#'

PDB_read <- function(PDB_ID){

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

}
