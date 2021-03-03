#' Perform basic PubMed search.
#'
#' @name pmtk_search_pubmed
#' 
#' @param pmed_search A character string, or list of character strings 
#' @param translate_syntax boolean: T to translate term to NCBI/mesh+tiab
#' @param verbose boolean: T to output progress
#' 
#' @return A data frame of PMIDs   
#' 
#' @export
#' @rdname pmtk_search_pubmed
#' 
pmtk_search_pubmed <- function (pmed_search,
                                translate_syntax = T,
                                verbose = T,
                                retmax = 5e6) {
  
  ## good eg: "violent depression"
  pre_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    
  returns <- lapply(1:length(pmed_search), function(y) {
    
    s1 <- pmed_search[y]
    
    if (!translate_syntax) {s2 <- s1} else{
      s2 <- paste0(s1, '[MH]', ' OR ', s1, '[TIAB]')
    }
    
    #nns <- rentrez::entrez_search(term = s2, db = 'pubmed')$count
    nns = format(retmax, scientific = F)
  
    if (nns == 0) {
      
      data.table::data.table(search = s1, pmid = NA)} else{ 
        
          url_term_query <- gsub(" ", "+", s1, fixed = TRUE)
          
          ## Make url
          full_url <- paste0 (pre_url, 
                              "db=pubmed", 
                              "&retmax=", nns, 
                              "&term=", 
                              url_term_query,
                              '%5BMH%5D+OR+',
                              url_term_query,
                              '%5BTIAB%5D')
                              #"&usehistory=n" 
          
          # x <- RCurl::getURL(full_url) 
          x <- httr::GET(full_url)
          x1 <- xml2::read_xml(x) 
          x2 <- xml2::xml_find_all(x1, './/Id')
          x3 <- xml2::xml_text(x2)
          
          if (length(x3) == 0) {x3 <- NA}
          
          out <- data.table::data.table(search = s1, pmid = x3) 
          
          if(verbose){
            finally =  print(paste0(y, ' / ', length(pmed_search), ' ',
                                    s2, ': ', nrow(out), ' records'))}
          
          out
      } 
    
    })
  
  data.table::rbindlist(returns)
}
