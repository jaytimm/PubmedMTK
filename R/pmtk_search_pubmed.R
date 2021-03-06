#' Perform basic PubMed search.
#'
#' @name pmtk_search_pubmed
#' @param search_term Query term as character string
#' @param fields PubMed fields to query
#' @param sleep In seconds.
#' @return A data frame of PMIDs   
#' 
#' @export
#' @rdname pmtk_search_pubmed
#' 
#' 
pmtk_search_pubmed <- function (search_term, 
                                fields = c('TIAB','MH'),
                                sleep = 1) { # max_n
  
  pre_url1 <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  pre_url2 <- "db=pubmed&retmax=5000000&term="
    
  url_term_query <- gsub(" ", "+", search_term, fixed = TRUE)
  
  if(is.null(fields)) { 
    fields3 <- url_term_query 
    out_search <- search_term} else{
      
      fields0 <- paste0('%5B', fields, '%5D')
      fields1 <- paste0(url_term_query, fields0)
      fields2 <- paste0(fields1, sep = '+OR+', collapse = '')
      fields3 <- gsub('\\+OR\\+$', '', fields2) 
    
      out_search <- paste0(search_term, 
                           paste0('[', fields, ']'), 
                           collapse = ' OR ')
    }
  
  ## Make url
  full_url <- paste0 (pre_url1, pre_url2, fields3)
  
  x <- httr::GET(full_url)
  x1 <- xml2::read_xml(x) 
  x2 <- xml2::xml_find_all(x1, './/Id')
  x3 <- xml2::xml_text(x2)
  
  if (length(x3) == 0) {x3 <- NA}
  
  out <- data.table::data.table(search_term = search_term, pmid = x3) 
  
  Sys.sleep(sleep)
  finally =  print(paste0(out_search, ': ', nrow(out), ' records'))
  out
}
