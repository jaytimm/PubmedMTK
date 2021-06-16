#' Search corpus for lexico-gram patterns in context via corpus::text_locate function.
#' Details.
#'
#' @name pmtk_locate_search
#' @param text character vector
#' @param doc_id character vector
#' @param search string
#' @param window integer
#' @param stem boolean
#' @return A data frame.
#'
#' @export
#' @rdname pmtk_locate_search
#'
pmtk_locate_search <- function(text,
                               doc_id,
                               search,
                               stem = F,
                               window = 15) {
  
  if(is.list(text)) {text <- unlist(lapply(text, paste0, collapse = ' '))}
  
  if(stem)(stx <- 'en') else{stx <- NULL}
  
  y <- corpus::text_locate(x = text, terms = search, stemmer = stx) 
  
  first2 <- sprintf("^\\S+( \\S+){0,%d}", window)

  y$rhs <- stringi::stri_extract(trimws(y$after), regex = first2)
  y$lhs <- stringi::stri_extract(stringi::stri_reverse(trimws(y$before)), regex = first2)
  y$lhs <- stringi::stri_reverse(y$lhs) ## -- ???!

  class(y) <- 'data.frame'
  y$text <- as.integer(y$text)
  y$doc_id <- doc_id[y$text]
  y$instance <- as.character(y$instance)
  
  y[, c('doc_id', 'lhs', 'instance', 'rhs')]
}
