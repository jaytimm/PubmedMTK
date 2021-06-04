#' Build data structures 2d visualizations.
#'
#' @name pmtk_2d
#' @param x A topic-term summary df.
#' @return A list. 
#' 
#' @export
#' @rdname pmtk_2d
#' 
pmtk_2d <- function(mat,
                    seed = 999,
                    d.method = "euclidean",
                    c.method = "ward.D") {  
  
  ##### sparse matrix
  # mat <- tidytext::cast_sparse(data = df,
  #                               row = topic_id,
  #                               column = feature, 
  #                               value = beta)
  
  ##### tSNE
  set.seed(seed)
  tsne <- Rtsne::Rtsne(X = as.matrix(mat), 
                       check_duplicates = T,
                       perplexity = 5)
  
  tsne1 <- data.frame(topic_id = rownames(mat), tsne$Y)
  
  ##### PCA
  cormat <- cor(data.frame(t(as.matrix(mat))))
  pca <- data.frame(prcomp(dist(cormat), 
                           scale = T, 
                           center = T)$x, 
                    row.names = NULL)
  
  pca1 <- pca[,c('PC1', 'PC2')]
  pca1$topic_id <- 1:ncol(cormat)
  colnames(pca1)[1:2] <- c('X1', 'X2')
  pca1 <- pca1[, c(3, 1:2)]
  
  ##### hierarchical cluster -- for dendrogram --
  disthi <- dist(mat, method = d.method) 
  hc1 <- hclust(disthi, method = c.method) 
  
  ##### output --
  list('pca' = pca1, 'tsne' = tsne1, 'hc' = hc1)
  
}
