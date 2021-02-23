#' Perform basic PubMed search.
#'
#' @name pmtk_get_pmids
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
#' @rdname pmtk_get_pmids
#' 
pmtk_get_pmids <- function (pmed_search,
                            search_as_is = T, ##
                            db = 'pubmed',  ##
                            max_url = 2e6,
                            summarize = F, ## - ??
                            verbose = F) {
  
  ## we need to clean up these parameters -- too many -- 
  
  x <- lapply(1:length(pmed_search), function(y) {
    
    s1 <- pmed_search[y]
    
    if (search_as_is) {s2 <- s1} else{
      s2 <- paste0(s1, '[MH]', ' OR ', s1, '[TIAB]')}
    # depression[MH] OR depression[tIAB]
    rentrez_search <- rentrez::entrez_search(term = s2, db = db)
    ## Needs to be 
    url_count <- min(max_url, rentrez_search$count) 
    
    if (verbose) {finally =  print(paste0(y, ' / ', length(pmed_search)))}
    
    if (url_count == 0) { 
      
      data.frame(search = s1, 
                 pmid = NA,
                 stringsAsFactors = FALSE)} else{
                   url_term_query <- gsub(" ", "+", s1, fixed = TRUE)
                   pmids <- paste0 ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
                                    "db=", db, "&retmax=", 
                                    url_count, 
                                    "&term=", 
                                    url_term_query, 
                                    "&usehistory=n") 
                   
                   x <- RCurl::getURL(pmids)
                   x1 <- xml2::read_xml(x) 
                   x2 <- xml2::xml_find_all(x1, './/Id')
                   x3 <- xml2::xml_text(x2)
                   
                   if (length(x3) == 0) {x3 <- NA}
                   data.frame(search = s1, pmid = x3) } })
  
  x <- data.table::rbindlist(x)
  
  if (summarize == FALSE) {x} else{
    x <- x[order(pmid, search)]
    x[, list(search = paste(search, collapse = " | ")), by = pmid]
  }
}
