#' Aggregate search results by query term combinations.
#'
#' @name pmtk_query_bigrams
#' @param search_results
#' @return A data frame  
#' 
#' @importFrom reshape2 melt
#' 
#' @export
#' @rdname pmtk_query_bigrams
#'
pmtk_query_bigrams <- function(search_results){
  
  search_results$value = 1 ### change such that we can add this as value-value -- 
  
  pmatrix <- tidytext::cast_sparse(data = search_results,
                                   row = pmid, 
                                   column = search, 
                                   value = value)
  
  v0 <- PubmedMTK::mtk_dtm_tcm(pmatrix)
  
  v1 <- cbind('search' = unique(search_results$search),
              as.data.frame(as.matrix(v0), 
                            row.names = NA))
  
  v2 <- reshape2::melt(v1, 'search', c(2:ncol(v1)))
  return(v2)
}
