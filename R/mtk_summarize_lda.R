#' Summarize topic model output
#' LDA model returned from `text2vec::LDA` - a fast/efficient LDA implementation - excellent for exploration -
#' 
#' @name mtk_summarize_lda
#' @param lda An lda object returned from `text2vec::LDA`
#' @param topic_feats_n Number of features to include in topic summary
#' @return 
#' 
#' @export
#' @rdname mtk_summarize_lda
#' 
mtk_summarize_lda <- function (lda, topic_feats_n = 10){
  
  ## Extract topic_word_distributions
  twd <- data.table::data.table(lda$.__enclos_env__$private$topic_word_distribution_with_prior())

  twd$topic_id <- 1:nrow(twd)
  twd1 <- data.table::melt.data.table(twd, id.vars = 'topic_id')
  data.table::setnames(twd1, old = "value", new = "beta")
  data.table::setnames(twd1, old = "variable", new = "feature")
  twd2 <- data.table::setorder(twd1,topic_id, -beta)
  twd3 <- twd2[, head(.SD, topic_feats_n), keyby = topic_id]
  tws <- twd3[ , .(topic_features = paste0(feature, collapse = ' | ')), by = topic_id]
  
  out <- list("topic_word_dist" = twd2, "topic_summary" = tws)
  out
}

  

  # ## Extract document_topic_distributions
  # dtds <- data.frame(lda$.__enclos_env__$private$doc_topic_distribution(), 
  #                    stringsAsFactors = FALSE) 
  # 
  # dtd <- cbind(doc_id = unique(dtm$doc_id), dtds)
  # dtd <- reshape2::melt(dtd, id.vars = 'doc_id')
  # dtd <- subset(dtd, value > 0)
  # dtd$variable <- gsub('X', '', dtd$variable)
  # 
  # ## setNames business here -- 
  # colnames(dtd)<- c('doc_id', 'topic', 'score')
  # 
  # 
  # ## labels
  # twd1 <- data.table::setorder(data.table::setDT(twd), -score)[, head(.SD, 4), keyby = topic]
  # twd1 <- twd1[ , .(topic_label = toString(descriptor_name)), by = topic]
  # # twd1$topic_label <- paste0(topic, ' - ', topic_label)
  # twd <- merge(twd, twd1, by = 'topic')
  # 
  # dtd <- merge(dtd, twd1, by = 'topic')
  # 
  # 
  # 
  # summary <- twd[ , .(topic_label = toString(descriptor_name)), 
  #                 by = topic]
  # 
  # out <- list("topic_word_dist" = twd, 
  #             "doc_topic_dist" = dtd,
  #             'topic_pca_coords' = topic_pca_coords,
  #             
  #             'lda_model' = lda, ## for viz -- ?? -- 
  #             'summary' = summary)
  # 
  # out }

