#' PDB Downloading and Reading
#'
#' @param PDB_ID PDB ID
#'
#' @return Matrix
#' @export
#'

PDB_read <- function(PDB_ID){

  pdb <- suppressMessages(bio3d::read.pdb(bio3d::get.pdb(PDB_ID, URL=TRUE)))

  atom_matrix <- pdb$atom
  atom_matrix <- atom_matrix[which(atom_matrix$type == "ATOM"),]
  atom_matrix <- atom_matrix[which(atom_matrix$resid != "HOH"),]

  message(crayon::white(paste0("STEP:","\n" ,"Reading PDB:", "\t\t\t", "DONE")))

  return(atom_matrix)
}
