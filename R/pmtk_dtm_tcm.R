#' Convet DTM to TCM; ie, treat document as context.
#'
#' @name pmtk_dtm_tcm
#' @param dtm 
#' @return A sparse matrix  
#' 
#' @importFrom Matrix t
#' 

#' @export
#' @rdname pmtk_dtm_tcm
#'
pmtk_dtm2tcm <- function(dtm){
  
  # create a binary matrix
  dtm_binary <- dtm > 0
  
  # dot product gives us the result
  result <- Matrix::t(dtm_binary) %*% dtm
  
  result
}