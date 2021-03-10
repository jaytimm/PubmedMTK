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
  
  ########
  
  x <- as(x, "character")
  newData <- XML::xmlParse(x)
  records <- XML::getNodeSet(newData, "//PubmedArticle")
  
  pmid <- XML::xpathSApply(newData,"//MedlineCitation/PMID", 
                           XML::xmlValue)
  
  year <- lapply(records, XML::xpathSApply, ".//PubDate/Year", 
                 XML::xmlValue) 
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
  
  Encoding(rownames(y)) <- 'UTF-8'    
  return(y)
}