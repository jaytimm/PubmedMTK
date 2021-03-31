#' Build data structures 2d visualizations.
#'
#' @name pmtk_2d
#' @param x A topic-term summary df.
#' @return A list. 
#' 
#' @export
#' @rdname pmtk_2d
#' 
pmtk_2d <- function(data,
                    row,
                    column, 
                    value, 
                    d.method = "euclidean",
                    c.method = "ward.D") {  
  
  ##### sparse matrix
  matx <- tidytext::cast_sparse(data = data,
                                row = row,
                                column = column, 
                                value = value)
  
  ##### tSNE
  set.seed(999)
  tsne <- Rtsne::Rtsne(X = as.matrix(matx), 
                       check_duplicates = FALSE,
                       perplexity = 5)
  
  tsne1 <- data.frame(topic_id = rownames(matx), tsne$Y)
  
  ##### PCA
  cormat <- cor(data.frame(t(as.matrix(matx))))
  pca <- data.frame(prcomp(dist(cormat), 
                           scale = T, 
                           center = T)$x, 
                    row.names = NULL)
  
  pca1 <- pca[,c('PC1', 'PC2')]
  pca1$topic_id <- 1:ncol(cormat)
  colnames(pca1)[1:2] <- c('X1', 'X2')
  pca1 <- pca1[, c(3, 1:2)]
  
  ##### hierarchical cluster -- for dendrogram --
  disthi <- dist(matx, method = d.method) 
  hc1 <- hclust(disthi, method = c.method) 
  
  ##### output --
  list('pca' = pca1, 'tsne' = tsne1, 'hc' = hc1)
  
}