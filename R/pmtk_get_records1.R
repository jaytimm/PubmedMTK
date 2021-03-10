#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_get_records1
#' @param x A character vector of PMIDs
#' @return A data frame.
#' 
#' 
#' @export
#' @rdname pmtk_get_records1
pmtk_get_records1 <- function (x) {

  fetch.pubmed <- rentrez::entrez_fetch(db = "pubmed", 
                                        id = x,
                                        rettype = "xml", 
                                        parsed = T)
  
  as(fetch.pubmed, "character")}
  ########
  
