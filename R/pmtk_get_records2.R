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
                             FUN = pmtk_get_records1,
                             cl = clust)
  
  parallel::stopCluster(clust)
  return(mess2)
}


#######
pmtk_get_records1 <- function (x) {

  x1 <- tryCatch({rentrez::entrez_fetch(db = "pubmed", 
                                        id = x,
                                        rettype = "xml",
                                        parsed = T) }, 
                 error = function(e) {'error-unspecified'} )
  
  Sys.sleep(0.25)
  if (any(!class(x1) == 'character')) {
    PubmedMTK:::pmtk_strip_xml(x1) } else{}
}



#######
pmtk_strip_xml <- function (x) {
  
  newData <- XML::xmlParse(as(x, "character"))
  records <- XML::getNodeSet(newData, "//PubmedArticle")
  
  pmid <- XML::xpathSApply(newData,"//MedlineCitation/PMID", XML::xmlValue)
  
  year <- lapply(records, XML::xpathSApply, ".//PubDate/Year", XML::xmlValue) 
  year[sapply(year, is.list)] <- NA
  year[which(sapply(year, is.na) == TRUE)] <- lapply(records[which(sapply(year, is.na) == TRUE)], 
                                                     XML::xpathSApply, ".//PubDate/MedlineDate", XML::xmlValue)
  
  ## as date formal -- !!! -- 
  year <- gsub(" .+", "", year)
  year <- gsub("-.+", "", year)
  
  articletitle <- lapply(records, XML::xpathSApply, ".//ArticleTitle", XML::xmlValue) 
  articletitle[sapply(articletitle, is.list)] <- NA  
  articletitle <- unlist(articletitle)
  
  meshHeadings <- lapply(records, XML::xpathSApply, ".//DescriptorName", XML::xmlValue)
  meshHeadings[sapply(meshHeadings, is.list)] <- NA
  meshHeadings <- sapply(meshHeadings, paste, collapse = "|")
  
  chemNames <- lapply(records, XML::xpathSApply, ".//NameOfSubstance", XML::xmlValue)
  chemNames[sapply(chemNames, is.list)] <- NA
  chemNames <- sapply(chemNames, paste, collapse = "|")
  
  keywords <- lapply(records, XML::xpathSApply, ".//Keyword", XML::xmlValue)
  keywords[sapply(keywords, is.list)] <- NA
  keywords <- sapply(keywords, paste, collapse = "|")
  
  y <- data.frame(pmid, year, articletitle,
                  meshHeadings, chemNames, 
                  keywords) 
  
  abstract <- lapply(records, XML::xpathSApply, ".//Abstract/AbstractText", XML::xmlValue)
  abstract[sapply(abstract, is.list)] <- NA
  abstract <- sapply(abstract, paste, collapse = " | ")
  y$abstract <- abstract   
  
  Encoding(rownames(y)) <- 'UTF-8'    
  
  data.table::setDT(y)
  clean_nas <- function(x) {
    ifelse(x %in% c(' ', 'NA', 'n/a', 'n/a.') | is.na(x), NA, x) }
  cols <- colnames(y)
  y[, c(cols) := lapply(.SD, clean_nas), .SDcols = cols]
  return(y)
}