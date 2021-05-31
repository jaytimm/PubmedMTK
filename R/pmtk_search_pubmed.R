#' Perform basic PubMed search.
#'
#' @name pmtk_search_pubmed
#' 
#' @param search_term A character string, or list of character strings 
#' @param translate_syntax boolean: T to translate term to NCBI/mesh+tiab
#' @param verbose boolean: T to output progress
#' 
#' @return A data frame of PMIDs   
#' 
#' @export
#' @rdname pmtk_search_pubmed
#' 
#' 

pmtk_search_pubmed <- function (search_term, 
                                fields = c('TIAB','MH'),
                                sleep = 1) { #
  
  ## good eg: "violent depression"
  pre_url1 <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  pre_url2 <- "db=pubmed&retmax=5000000&term="
    
  #s2 <- paste0(search, '[MH]', ' OR ', s1, '[TIAB]')
  #nns <- rentrez::entrez_search(term = s2, db = 'pubmed')$count
  #nns <- format(5e6, scientific = F)

  url_term_query <- gsub(" ", "+", search_term, fixed = TRUE)
  
  if(is.null(fields)) { fields3 <- url_term_query } else{
    fields0 <- paste0('%5B', fields, '%5D')
    fields1 <- paste0(url_term_query, fields0)
    fields2 <- paste0(fields1, sep = '+OR+', collapse = '')
    fields3 <- gsub('\\+OR\\+$', '', fields2) }
  
  ## Make url
  full_url <- paste0 (pre_url1, pre_url2, fields3)
  
  x <- httr::GET(full_url)
  x1 <- xml2::read_xml(x) 
  x2 <- xml2::xml_find_all(x1, './/Id')
  x3 <- xml2::xml_text(x2)
  
  if (length(x3) == 0) {x3 <- NA}
  
  out <- data.table::data.table(search_term = search_term, pmid = x3) 
  
  Sys.sleep(sleep)
  finally =  print(paste0(search_term, ': ', nrow(out), ' records'))
  out
}
