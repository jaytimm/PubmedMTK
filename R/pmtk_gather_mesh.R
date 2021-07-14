#' Extract MeSH terms.
#'
#' @name pmtk_gather_mesh
#' @param x A data frame returned from pmtk_get_records2()
#' @return A data frame 
#' 
#' @export
#' @rdname pmtk_gather_mesh
#' 
pmtk_gather_mesh <- function (x) {

  data.table::setDT(x)
  x0 <- data.table::melt(x, 
                         measure.vars = c('keywords', 
                                          'meshHeadings',
                                          'chemNames'),
                         variable.name = "type", 
                         value.name = "term",
                         na.rm = T)
  
  x0[, term := tolower(trimws(term))]
  x1 <- x0[, list(term = unlist(strsplit(term, "\\|"))), 
           by = list(pmid, type)]
  x1[, term := gsub(' ', '_', term)]
  
  x1[, type := factor(type, levels = c('meshHeadings', 'chemNames', 'keywords'))]
  data.table::setorder(x1, pmid, type)
  
  unique(x1, by = c('pmid', 'term'))
}
