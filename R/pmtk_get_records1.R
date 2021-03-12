#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_get_records1
#' @param x A character vector of PMIDs
#' @return A data frame
#' 
#' 
#' @export
#' @rdname pmtk_get_records1
pmtk_get_records1 <- function (x) {
  
  ## useful error message here -- 
  
  if(length(x) > 199) {
    message('Max 199 records per query -- only first 199 returned')
    x <- x[1:199] 
    }
  
  x1 <- rentrez::entrez_fetch(db = "pubmed", 
                              id = x,
                              rettype = "xml", 
                              parsed = T)
  
  pmtk_strip_xml1(x1)
  }
  


## https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html


## Information found in PubMed that indicates it is "indexed by MEDLINE" is considered peer reviewed. Look for the phrase "indexed by MEDLINE" under the citation or abstract information.

## PubMed does not provide a search filter to limit to only peer reviewed articles. For other citations, look up the journal title in the NCBI Journals Database, click on the journal title, find a publisher's website link and go to that website. Look for something on the page that gives details about the journal and then read through it to find if the journal goes through a peer review process. 


pmtk_strip_xml1 <- function (x,
                             include_abstract = F) {
  
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
  
  if(include_abstract) {
    
    abstract <- lapply(records, XML::xpathSApply, ".//Abstract/AbstractText", XML::xmlValue)
    abstract[sapply(abstract, is.list)] <- NA
    
    y$abstract <- abstract  }
  
  Encoding(rownames(y)) <- 'UTF-8'    
  return(y)
}
