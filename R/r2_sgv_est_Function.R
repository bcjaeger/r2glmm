

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
#' @param blksizes vector of block sizes
#' @param vmat The blocked covariance matrix
#' @return The SGV of the covariance matrix \code{vmat}.
#' @examples
#' library(Matrix)
#' v1 = matrix(c(1,0.5,0.5,1), nrow = 2)
#' v2 = matrix(c(1,0.2,0.1,0.2,1,0.3,0.1,0.3,1), nrow = 3)
#' v3 = matrix(c(1,0.1,0.1,0.1,1,0.2,0.1,0.2,1), nrow = 3)
#' calc_sgv(nblocks = 3, blksizes = c(2,3,3), vmat = Matrix::bdiag(v1,v2,v3))
#' @export calc_sgv

calc_sgv <- function(nblocks=NULL, blksizes = NULL, vmat){

  # lme allows the user to input a list of matrices directly
  # if this is done, computation-time is considerably lessened

  if(!inherits(vmat, 'list')){

    # initialize a matrix list and a starting point

    mlist = list()
    start = 1


    if(is.null(blksizes)){

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

      for(i in 1:nblocks){

        # Get a 'block' out of the covariance matrix

        start = ifelse(i==1, 1, 1+sum(blksizes[1:i-1]))
        end   = ifelse(i==nblocks, nrow(vmat), start-1+blksizes[i])

        # Pick out the i^th block
        mlist[[i]] = as.matrix(vmat[start:end, start:end])

      }

    }

  } else {

    mlist = vmat

  }

  sgv = try(lapply(mlist, function(mat){ log(det(mat)) / nrow(mat)}))

  if(inherits(sgv, 'try-error')){

    stop('SGV is non-finite. Consider rescaling outcomes and Predictors')

  }

  # Back transform from log scale after taking mean of log sgvs
  # This gives the geometric mean of the sgv of the blocks
  # One numerical complication is non-finite sgv values.

  sgvec = unlist(sgv)
  keep = is.finite(sgvec)

  if(length(sgvec) - sum(keep) > 0){
    warning("Some SGV estimates are non-finite and have been adjusted")

    sgvec[!keep] = max(sgvec[keep])

  }


  SGV = exp( mean(sgvec) )

  return(as.numeric(SGV))

}


#' pqlmer
#' @description Fit a GLMM model with multivariate normal random effects using Penalized Quasi-Likelihood for mermod objects.
#'
#' @param formula The lme4 model formula.
#' @param family a family function of the error distribution and link function to be used in the model.
#' @param data the dataframe containing the variables in the model.
#' @param niter Maximum number of iterations to perform.
#' @param verbose if TRUE, iterations are printed to console.
#' @return A pseudo linear mixed model of class "lme" .
#' @examples
#' # Compare lmer PQL with lme PQL
#'
#' library(MASS)
#'
#' lmePQL = glmmPQL(y ~ trt + week + I(week > 2), random = ~ 1 | ID,
#'                   family = binomial, data = bacteria,
#'                   verbose = FALSE)
#'
#' merPQL= pqlmer(y ~ trt + week + I(week > 2) + (1 | ID),
#'                family = binomial, data = bacteria,
#'                verbose = FALSE)
#'
#' summary(lmePQL)
#' summary(merPQL)
#' @export pqlmer

pqlmer <- function(formula, family, data, niter = 40, verbose = T){

    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("family' not recognized")
    }

  fit0 = stats::glm(lme4::nobars(formula), family = family, data = data)

  # Remove missing values (glm does this automatically)
  # data = data[!is.na(data[[as.character(fit0$formula[[2]])]]),]

  eta <- fit0$linear.predictors
  wz=NULL

  data$zz <- eta + fit0$residuals
  data$wz <- fit0$weights

  fail_message = paste('Model failed to converge in', niter, 'iterations.')

  for (i in seq_len(niter)) {

    if (verbose) message(gettextf("iteration %d", i), domain = NA)

    if(i == niter) stop(fail_message)

    fit <- lme4::lmer(formula = stats::update(formula, zz ~ .),
                      data = data, weights = wz)

    etaold <- eta
    eta    <- stats::fitted(fit)

    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) break

    # mu = h(XB + Zb), where h is inverse link function
    mu <- family$linkinv(eta)

    # delta = derivative of mu with respect to eta
    mu.eta.val <- family$mu.eta(eta)

    data$zz <- eta + (fit0$y - mu)/mu.eta.val

    data$wz <- mu.eta.val^2/family$variance(mu)

  }

  return(fit)

}



