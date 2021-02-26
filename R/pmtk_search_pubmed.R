#' Perform basic PubMed search.
#'
#' @name pmtk_search_pubmed
#' 
#' @param pmids A vector of PMIDs 
#' @param out_file File path for output
#' @param file_prefix String specifying how batch file name
#' 
#' @return A data frame of PMIDs   
#' 
#' @importFrom rentrez entrez_search
#' @importFrom RCurl getURL
#' @import data.table xml2 
#' 
#' 
#' @export
#' @rdname pmtk_search_pubmed
#' 
pmtk_search_pubmed <- function (pmed_search,
                                translate_syntax = T,
                                verbose = T) {
  
  ## we need to clean up these parameters -- too many -- 
  db <- 'pubmed'
  #convert_syntax = T
  pre_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    
  returns <- lapply(1:length(pmed_search), function(y) {
    
    s1 <- pmed_search[y]
    
    if (!translate_syntax) {s2 <- s1} else{
      s2 <- paste0(s1, '[MH]', ' OR ', s1, '[TIAB]')
      }
    # depression[MH] OR depression[tIAB]
    rentrez_search <- rentrez::entrez_search(term = s2, db = db)
    
    # url_count <- min(max_url, rentrez_search$count) 
    url_count <- rentrez_search$count 
    
    if(verbose){
      finally =  print(paste0(y, ' / ', length(pmed_search), ' ',
                              s2, ': ', url_count, ' records'))}
    
    if (url_count == 0) { 
      
      data.table::data.table(search = s1, pmid = NA)} else{ 
                   
        ## below is different than rentrez_search
          url_term_query <- gsub(" ", "+", s1, fixed = TRUE)
          
          pmids <- paste0 (pre_url, 
                           "db=", db, "&retmax=", 
                           url_count, 
                           "&term=", 
                           url_term_query,
                           '%5BMH%5D+OR+',
                           url_term_query,
                           '%5BTIAB%5D')
                           #"&usehistory=n" 
          
          x <- RCurl::getURL(pmids)
          
          ## for search = induced senescence --??
          x1 <- xml2::read_xml(x) 
          x2 <- xml2::xml_find_all(x1, './/Id')
          x3 <- xml2::xml_text(x2)
          
          if (length(x3) == 0) {x3 <- NA}
          
          data.table::data.table(search = s1, pmid = x3) 
      
      } 
    
    })
  
  data.table::rbindlist(returns)
}
