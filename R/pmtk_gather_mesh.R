#' Download abstract and meta data for research articles included in PubMed.
#'
#' @name pmtk_gather_mesh
#' @param meta_df A vector of PMIDs 
#' @return A data frame 
#' 
#' @import data.table
#' 

## Extract KEYWORDS, MeSH HEADINGS & CHEM-NAMES from columns included in metadata as a clean data table. 

#' @export
#' @rdname pmtk_gather_mesh
#' 
pmtk_gather_mesh <- function (meta_df) {
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
  terms <- subset(terms, !term %in% c('', 'na'))
  out <- data.table::setDT(terms)[, list(term = unlist(strsplit(term, "\\|"))), 
                                  by = list(pmid, type)]

  out <- out[, list(cooc =.N), by = list(pmid, type, term)]
  colnames(out) <- c('pmid', 'type', 'descriptor_name', 'count')
  out
}