#' Load output from `pmtk_batch_abstracts` as a single data.frame
#' 
#' @name pmtk_strip_batches
#' @param in_file File path to batches
#' @param file_prefix  String specifying batch name
#' @return A data.frame  
#' @importFrom data.table rbindlist
#' 
#' 
#' @export
#' @rdname pmtk_strip_batches
#' 
pubmed_strip_batches <- function (in_file,
                                 file_prefix) {
  setwd(in_file)
  gfiles <- list.files(path = in_file, 
                       pattern = file_prefix, 
                       recursive=TRUE) 
  ref <- lapply(gfiles, readRDS) 
  ref <- data.table::rbindlist(ref)
  
  tif <- ref[, c('pmid', 'text')]  #}
  
  meta <- ref[, c(setdiff(colnames(ref), colnames(tif)[-1])), 
              with = F]  
  
  list("meta" = meta, "tif" = tif)
}
