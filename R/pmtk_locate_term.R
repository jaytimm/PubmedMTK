#' Search corpus for lexical patterns in context via corpus::text_locate function.
#'
#' @name pmtk_locate_term
#' @param text character vector
#' @param doc_id character vector
#' @param term string
#' @param window integer
#' @param stem boolean
#' @return A data frame.
#'
#' @export
#' @rdname pmtk_locate_term
#'
pmtk_locate_term <- function(text,
                               doc_id,
                               term,
                               stem = F,
                               window = 15) {
  
  if(is.list(text)) {text <- unlist(lapply(text, paste0, collapse = ' '))}
  
  if(stem)(stx <- 'en') else{stx <- NULL}
  
  y <- corpus::text_locate(x = text, 
                           terms = term, 
                           stemmer = stx) 
  
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
