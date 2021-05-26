#' Structure PMC article as data frame, grouped by article section
#'
#' @name pmtk_extract_pmc
#' @param x A file path to a PMC article
#' @return A data frame 
#' 
#' @export
#' @rdname pmtk_extract_pmc
#' 
# ldir <- '/home/jtimm/Desktop/oa_bulk/xml-demo/'
# files <- dir(ldir, full.names = TRUE, recursive = T)

pmtk_extract_pmc <- function(x) {
  x0 <- xml2::read_xml(x)
  x1 <- xml2::xml_child(x0, 2)
  
  header_titles <- lapply(xml2::xml_children(x1),
                          function(x) {
                            xml2::xml_text(xml2::xml_find_first(x, ".//title"))}
  )
  
  sub_titles <- lapply(xml2::xml_children(x1),
                       function(x) {
                         xml2::xml_text(xml2::xml_find_all(x, ".//title"))}
  )
  
  toc <- setNames(stack(setNames(sub_titles, header_titles)), c('sub', 'header'))
  toc0 <- subset(toc, sub != header)
  
  text <- lapply(xml2::xml_children(x1), xml2::xml_text)
  
  df <- data.frame(pmid = gsub('(^.*)(PMC)([0-9]*)(\\.nxml$)', '\\2\\3', x),
                   section = unlist(header_titles), 
                   text = unlist(text),
                   row.names = NULL)
  
  heads <- paste0('^', df$section, collapse = '|')
  heads <- gsub('\\(', '\\\\(', heads)
  heads <- gsub('\\)', '\\\\)', heads)
  
  df$text <- gsub(heads, '', df$text)
  
  if (nrow(toc0) > 0) {
    splits <- paste0(toc$sub, collapse = '|')
    splits <- gsub('\\(', '\\\\(', splits)
    splits <- gsub('\\)', '\\\\)', splits)
    q <- paste0('(', splits, ')')
    df$text <- gsub(q, '\n\n\\1: ', df$text)
  }
  
  df
}

## xx <- lapply(files[10:20], pmtk_extract_pmc)