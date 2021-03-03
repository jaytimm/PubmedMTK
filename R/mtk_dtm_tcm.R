#' Convet DTM to TCM; ie, treat document as context.
#'
#' @name mtk_dtm_tcm
#' @param dtm A document-term matrix
#' @return A sparse matrix  
#' 

#' @export
#' @rdname mtk_dtm_tcm
#'
mtk_dtm_tcm <- function(dtm){
  
  # create a binary matrix
  dtm_binary <- dtm > 0
  
  # dot product gives us the result
  result <- Matrix::t(dtm_binary) %*% dtm
  
  result
}