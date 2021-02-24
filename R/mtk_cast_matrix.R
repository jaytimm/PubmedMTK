#' Project data frame as a sparse matrix.
#'
#' @name mtk_cast_matrix
#' @param dtm 
#' @param row
#' @param column
#' @param value
#' @return A sparse matrix  
#' 
#' @importFrom Matrix sparseMatrix
#' 

#' @export
#' @rdname mtk_cast_spmatrix
#'
mtk_cast_matrix <- function(dtm, row, column, value, ...) {
  
  ## -- ??
  row_col <- quo_name(enquo(row))
  column_col <- quo_name(enquo(column))
  value_col <- enquo(value)
  
  
  row_names <- data[[row_col]]
  col_names <- data[[column_col]]
  if (is.numeric(value_col)) {
    values <- value_col
  } else {
    value_col <- quo_name(value_col)
    values <- data[[value_col]]
  }
  
  row_u <- unique(row_names)
  i <- match(row_names, row_u)
  col_u <- unique(col_names)
  j <- match(col_names, col_u)
  
  ret <- Matrix::sparseMatrix(
    i = i, j = j, x = values,
    dimnames = list(row_u, col_u), ...
  )
  
  ret
}