

#' Checks if a matrix is Compound Symmetric.
#'
#' @param mat The matrix to be tested.
#' @param tol a number indicating the smallest acceptable difference between off diagonal values.
#' @return True if the matrix is compound symmetric.
#' @examples
#'
#' gcmat <- matrix(c(1,0.2,0.1,0.2,1,0.3,0.1,0.3,1), nrow = 3)
#' csmat <- matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow = 3)
#' is.CompSym(csmat)
#' @export is.CompSym

is.CompSym = function(mat, tol = 0.00001){

  if(length(mat)==1) {

    return(FALSE)

  } else {

    off.diag.vals = mat[upper.tri(mat)]
    isTRUE(abs(max(off.diag.vals) - min(off.diag.vals)) < tol)

  }

}

#' Compute the standardized generalized variance (SGV) of a blocked diagonal matrix.
#'
#' @param nblocks Number of blocks in the matrix.
#' @param vmat The blocked covariance matrix
#' @return The SGV of the covariance matrix \code{vmat}.
#' @examples
#' library(Matrix)
#' v1 = matrix(c(1,0.5,0.5,1), nrow = 2)
#' v2 = matrix(c(1,0.2,0.1,0.2,1,0.3,0.1,0.3,1), nrow = 3)
#' calc_sgv(nblocks = 2, vmat = Matrix::bdiag(v1,v2))
#' @export calc_sgv

calc_sgv <- function(nblocks=NULL, vmat){

  # lme allows the user to input a list of matrices directly
  # if this is done, computation-time is considerably lessened

  if(class(vmat) != 'list'){

    # initialize a matrix list and a starting point

    mlist = list()
    start = 1

    # loop through the matrix blocks

    for(i in seq(nblocks)){

      # the stopping point is determined by the first zero
      # that follows after the positive real numbers
      # in the variance-covariance matrix. At each iteration,
      # update the starting and stopping points using previous

      stop = ifelse(i < nblocks,
                    which(vmat[start, start:ncol(vmat)]==0)[1] + (start-2),
                    nrow(vmat))

      mlist[[i]] = as.matrix(vmat[start:stop,start:stop])

      start = stop + 1

    }

  } else {

    mlist = vmat

  }

  sgv = try(lapply(mlist, function(mat){ log(det(mat)) / nrow(mat)}))

  if(class(sgv)=="try-error"){

    stop('SGV is non-finite. Consider rescaling outcomes and Predictors')

  }

  # Back transform from log scale after taking mean of log sgvs
  # This gives the geometric mean of the sgv of the blocks

  SGV = exp(mean(unlist(sgv)))

  return(as.numeric(SGV))

}

#' r2.sgv
#'
#' @description  Computes R squared for a linear or generalized linear mixed
#' model using penalized quasi likelihood (PQL) estimation and standardized
#' generalized variance (SGV). Currently implemented for linear mixed models
#' with \code{\link{lmer}} and \code{\link{lme}} objects.
#' For generalized linear mixed models, only \code{\link{glmmPQL}}
#' objects are compatible.
#'
#' @param model a fitted mixed model.
#' @param adj if TRUE, an adjusted R squared for model selection is provided.
#' @return The R squared statistic and adjusted R squared based on SGV
#'         for \code{model}.
#' @examples
#'
#' library(nlme)
#' m = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#' r2.sgv(m)
#'
#' library(MASS)
#' PQL_mod = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#' family = binomial, data = bacteria, verbose = F)
#' r2.sgv(PQL_mod)







