#' A sentence tokenizer that mostly works.
#'
#' @name pmtk_toke_sentences
#' @param text character vector
#' @param doc_id character vector
#' @return A data frame.
#'
#' @export
#' @rdname pmtk_toke_sentences
#'
pmtk_toke_sentences <- function(text,
                                doc_id){

  y <- text
  names(y) <- doc_id

  corp0 <- corpus::text_split(y,
                              units = 'sentences',
                              filter = corpus::text_filter(sent_suppress = c(corpus::abbreviations_en)))

  class(corp0) <- 'data.frame'
  corp0$text <- as.character(corp0$text)
  corp0$doc_id <- paste0(corp0$parent, '.', corp0$index)
  corp1 <- corp0[, c(4, 3)]
}

## xx <- sentence_tokenizer(x = corpus1)