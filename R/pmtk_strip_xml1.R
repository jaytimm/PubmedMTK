#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_strip_xml1
#' @param x XML output from NCBI API query. 
#' @return A list of data frames
#' 
#' @export
#' @rdname pmtk_strip_xml1
#' 
pmtk_strip_xml1 <- function (x) {
  
  newData <- XML::xmlParse(z)
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
  
  #Sys.sleep(0.1)
  Encoding(rownames(y)) <- 'UTF-8'    
  return(y)
}
