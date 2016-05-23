
#' Generate partial contrast matrices
#'
#' @param rows Number of rows in the contrast matrix
#' @param cols Number of columns in the contrast matrix
#' @param index A number corresponding to the position of the fixed effect
#'  in the vector of fixed effect parameter estimates.
#' @return A contrast matrix designed to test the fixed effect
#'  corresponding to \code{index} in the vector of fixed effects.
#' @examples
#' make.partial.C(4, 5, 2)
#' make.partial.C(4, 5, 3)
#' @export make.partial.C
make.partial.C = function(rows, cols, index){
  x=matrix(0, nrow = rows, ncol = cols)
  x[index-1, index] = 1
  return(x)
}
