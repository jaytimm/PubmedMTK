#' Load output from `pmtk_download_abs` as a single data.frame
#' 
#' @name pmtk_loadr_abs
#' @param in_file File path to batches
#' @param file_prefix  String specifying batch name
#' @return A data.frame  
#' 
#' @export
#' @rdname pmtk_loadr_abs
#' 
pmtk_loadr_abs <- function (in_file,
                            file_prefix) {

  setwd(in_file)
  gfiles <- list.files(path = in_file, 
                       pattern = file_prefix, 
                       recursive=TRUE) 
  ref <- lapply(gfiles, readRDS) 
  ref <- data.table::rbindlist(ref)
  
  clean_nas <- function(x) {
    ifelse(x %in% c(' ', 'NA', 'n/a', 'n/a.') | is.na(x), NA, x) }
  
  cols <- colnames(ref)
  ref[, c(cols) := lapply(.SD, clean_nas), .SDcols = cols]
  
  tif <- ref[, c('pmid', 'text')]  #}
  
  meta <- ref[, c(setdiff(colnames(ref), colnames(tif)[-1])), with = F]  
  
  list("meta" = meta, "tif" = tif)
}


