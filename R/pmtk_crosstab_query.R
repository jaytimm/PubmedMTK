#' Aggregate search results by query term combinations.
#'
#' @name pmtk_crosstab_query
#' @param search_results
#' @return A data frame  
#' 
#' @importFrom reshape2 melt
#' 
#' @export
#' @rdname pmtk_crosstab_query
#'
pmtk_crosstab_query <- function(search_results){
  
  search_results$value = 1 
  pmatrix <- tidytext::cast_sparse(data = search_results,
                                   row = pmid, 
                                   column = search, 
                                   value = value)
  
  v0 <- PubmedMTK::mtk_dtm_tcm(pmatrix)
  
  v1 <- cbind('search' = unique(search_results$search),
              as.data.frame(as.matrix(v0), 
                            row.names = NA))
  
  v2 <- reshape2::melt(v1, 'search', c(2:ncol(v1)))
  colnames(v2) <- c('term1', 'term2', 'n1n2')
  
  v2
  ## 
  v3 <- subset(v2, term1 == term2)
  v4 <- merge(v2, data.frame(term1 = v3$term1, n1 = v3$n))
  v5 <- merge(v4, data.frame(term2 = v3$term1, n2 = v3$n))
  
  v5$term2 <- as.character(v5$term2)
  v5 <- v5[order(v5$term1),]
  v5 <- subset(v5, term1 != term2)
  v5[, c(2, 1, 4, 5, 3)]
}
