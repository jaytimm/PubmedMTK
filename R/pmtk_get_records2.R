#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_get_records2
#' @param pmids A vector of PMIDs 
#' @param cores Numeric specifying number of cores to use 
#' @return A list of data frames
#' 
#' @export
#' @rdname pmtk_get_records2
#' 
#' 
pmtk_get_records2 <- function (pmids, 
                               cores, 
                               ncbi_key = NULL) {
  
  if(is.null(ncbi_key) & cores > 3) cores <- min(parallel::detectCores() - 1, 3)
  if(!is.null(ncbi_key)) rentrez::set_entrez_key(ncbi_key)
  
  batches <- split(pmids, ceiling(seq_along(pmids)/199)) 
  
  clust <- parallel::makeCluster(cores)
  #parallel::clusterEvalQ(cl = clust, library(PubmedMTK))
  parallel::clusterExport(cl = clust, 
                          varlist = c('batches'))
  
  
  mess2 <- pbapply::pblapply(X = batches,
                             FUN = pmtk_get_records1,
                             cl = clust)
  
  parallel::stopCluster(clust)
  return(mess2)
}
