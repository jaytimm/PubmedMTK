#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_get_records
#' @param pmids A vector of PMIDs 
#' @param cores Numeric specifying number of cores to use 
#' @return A list of data frames
#' 
#' @export
#' @rdname pmtk_get_records
#' 
#' 
pmtk_get_records <- function (pmids, cores) {
  
  batches <- split(pmids, ceiling(seq_along(pmids)/199)) 
  
  clust <- parallel::makeCluster(cores)
  parallel::clusterExport(cl = clust, 
                          varlist = c('x_get_records', 
                                      'x_strip_xml')) ## -- ??
  mess2 <- pbapply::pblapply(X = batches,
                             FUN = x_get_records,
                             cl = clust)
  parallel::stopCluster(clust)
  return(mess2)
}



x_get_records <- function (x) {
  
  fetch.pubmed <- rentrez::entrez_fetch(db = "pubmed", 
                                        id = x,
                                        rettype = "xml", 
                                        parsed = T)
  
  x1 <- as(fetch.pubmed, "character")
  x2 <- x_strip_xml(x1)
  Encoding(rownames(x2)) <- 'UTF-8'    
  return(x2)
}


x_strip_xml <- function(x) {
  
  ## journal, is_peer_reviewed + param -- include_abstract -- ??
  
  
  newData <- XML::xmlParse(x)
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

  data.frame(pmid, year, articletitle,
             meshHeadings, chemNames, 
             keywords)
}

