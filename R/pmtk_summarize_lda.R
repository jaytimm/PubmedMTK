#' Summarize topic model output returned from text2vec::LDA()
#' 
#' @name pmtk_summarize_lda
#' @param lda An lda object returned from text2vec::LDA()
#' @param topic_feats_n Number of features to include in topic summary
#' @return A list of data frames
#' 
#' @export
#' @rdname pmtk_summarize_lda
#' 
pmtk_summarize_lda <- function (lda, topic_feats_n = 10){
  
  ## Extract document_topic_distributions
  dtds <- data.frame(lda$.__enclos_env__$private$doc_topic_distribution())
  dtd <- cbind(doc_id = names(lda$.__enclos_env__$private$doc_len), dtds)
  
  data.table::setDT(dtd)
  dtd <- data.table::melt.data.table(dtd, 
                                     id.vars = 'doc_id', 
                                     variable.name = 'topic', 
                                     value.name = 'beta')
  dtd <- subset(dtd, beta > 0)
  dtd$topic <- gsub('X', '', dtd$topic)  
  
  
  ## Extract topic_word_distributions
  twd <- data.table::data.table(lda$.__enclos_env__$private$topic_word_distribution_with_prior())
  twd$topic_id <- 1:nrow(twd)
  twd1 <- data.table::melt.data.table(twd, id.vars = 'topic_id', 
                                      variable.name = 'feature', 
                                      value.name = 'beta')
  twd2 <- data.table::setorder(twd1, topic_id, -beta)
  
  twd3 <- twd2[, head(.SD, topic_feats_n), keyby = topic_id]
  tws <- twd3[ , .(topic_features = paste0(feature, collapse = ' | ')), by = topic_id]
  
  list("topic_word_dist" = twd2, 
       "topic_summary" = tws,
       "doc_topic_dist" = dtd)
}
