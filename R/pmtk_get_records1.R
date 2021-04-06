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
  
  PubmedMTK:::pmtk_strip_xml(x1)
  }
  


## https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html


## Information found in PubMed that indicates it is "indexed by MEDLINE" is considered peer reviewed. Look for the phrase "indexed by MEDLINE" under the citation or abstract information.

## PubMed does not provide a search filter to limit to only peer reviewed articles. For other citations, look up the journal title in the NCBI Journals Database, click on the journal title, find a publisher's website link and go to that website. Look for something on the page that gives details about the journal and then read through it to find if the journal goes through a peer review process. 
