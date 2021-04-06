#' Extract MeSH terms.
#'
#' @name pmtk_gather_meshb
#' @param meta_df Metadata data frame returned from `pmtk_loadr_abs'` 
#' @return A data frame 
#' 
#' @export
#' @rdname pmtk_gather_meshb
#' 
pmtk_gather_meshb <- function (meta_df) {
  # specific to `pubmed_strip_batches` output --
  
  meta_df$keywords <- tolower(trimws(meta_df$keywords))
  meta_df$meshHeadings <- tolower(trimws(meta_df$meshHeadings))
  meta_df$chemNames <- tolower(trimws(meta_df$chemNames))
  
  terms <- rbind(meta_df[, c('pmid', 'keywords')],
                 meta_df[, c('pmid', 'meshHeadings')],
                 meta_df[, c('pmid', 'chemNames')],
                 use.names=FALSE)
  
  terms$type <- c(rep('keyword', nrow(meta_df)),
                  rep('mesh_heading', nrow(meta_df)),
                  rep('chem_name', nrow(meta_df)))
  
  colnames(terms)[2] <- 'term'
  
  ## -- ??
  #terms <- subset(terms, !term %in% c('', 'na'))
  
  
  out <- data.table::setDT(terms)[, list(term = unlist(strsplit(term, "\\|"))), 
                                  by = list(pmid, type)]

  out <- out[, list(cooc =.N), by = list(pmid, type, term)]
  # colnames(out) <- c('pmid', 'type', 'descriptor_name', 'count')
  # out

  
  ### need to resolve the TYPE ISSUE -- !! we did somewhere -- 
  tidytext::cast_sparse(data = out,
                        row = pmid,
                        column = term, 
                        value = cooc)
}
