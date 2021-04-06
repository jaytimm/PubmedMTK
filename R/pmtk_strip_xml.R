#' Extract text from xml nodes.
#' 
#' @name pmtk_strip_xml
#' @param x A character vector of PMIDs
#' @return A data frame

#' @export
#' @rdname pmtk_strip_xml
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