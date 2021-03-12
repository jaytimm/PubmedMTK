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
                               cores = 3, 
                               ncbi_key = NULL) {
  
  if(is.null(ncbi_key) & cores > 3) cores <- min(parallel::detectCores() - 1, 3)
  if(!is.null(ncbi_key)) rentrez::set_entrez_key(ncbi_key)
  
  batches <- split(pmids, ceiling(seq_along(pmids)/199)) 
  
  clust <- parallel::makeCluster(cores)
  parallel::clusterExport(cl = clust, 
                          varlist = c('batches'),
                          envir = environment())
  
  mess2 <- pbapply::pblapply(X = batches,
                             FUN = PubmedMTK::pmtk_get_records1,
                             cl = clust)
  
  parallel::stopCluster(clust)
  return(mess2)
}