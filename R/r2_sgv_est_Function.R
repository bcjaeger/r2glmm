



#' Compute the standardized generalized variance (SGV) of a blocked diagonal matrix.
#'
#' @param nblocks Number of blocks in the matrix.
#' @param blk.sizes  A vector containing the dimension of each block.
#' @param vmat The blocked covariance matrix
#' @return The SGV of the covariance matrix \code{vmat}.
#' @examples
#' v1 = matrix(c(1,0.5,0.5,1), nrow = 2)
#' v2 = matrix(c(1,0.2,0.1,0.2,1,0.3,0.1,0.3,1), nrow = 3)
#' calc.sgv(nblocks = 2, blk.sizes = c(2,3), vmat = bdiag(v1,v2))
#' @export calc.sgv

calc.sgv <- function(nblocks, blk.sizes, vmat){

  if(!require("Matrix",quietly=T)) stop("package 'Matrix' is essential")

  # initialize the sgv values
  sgv = rep(0, nblocks)

  for(i in 1:nblocks){

    # Get the first 'block' out of the covariance matrix
    if(i == 1){
      start = 1; end = blk.sizes[1]
      mat = vmat[start:end, start:end]
    }

    # if i >= 2, pick out the i^{th} block from vmat
    if(i >= 2){
      start = 1+sum(blk.sizes[1:i-1])
      end = start-1+blk.sizes[i]
      mat = vmat[start:end, start:end]
    }
    # determinants can be large, so log function is used for numeric stability
    mat = as.matrix(mat)
    sgv[i] = log(det(mat))/nrow(mat)
  }

  # Back transform from log scale after taking mean of log sgvs
  # This gives the geometric mean of the sgv of the blocks
  SGV = exp(mean(sgv))

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
#' @export

r2.sgv <- function(model, adj = T){

  if(!require("Matrix",quietly=T)) stop("package 'Matrix' is essential")
  if(!require("mgcv",quietly=T)) stop("package 'mgcv' is essential")
  if(!require("AICcmodavg",quietly=T)) stop("package 'AICcmodavg' is essential")

  # Get fixed effects
  beta = fixef(model)
  p <- length(beta)

  if(any(c('merModLmerTest', 'lmerMod') %in% class(model))){

    s2e = getME(model, 'sigma')^2
    X = getME(model, 'X')
    n = nrow(X)
    Z = as.matrix(lme4::getME(model, 'Z'))
    G = as.matrix(s2e * getME(model, 'Lambda') %*% getME(model, 'Lambdat'))

    # Covariance Matrix from the model
    SigHat = Z%*%G%*%t(Z) + s2e*diag(nrow(Z))

    clust.id = names(model@flist)[ length(model@flist) ]
    obsperclust = as.numeric(table(model@frame[,clust.id]))
    nclusts = length(obsperclust)

    SGV = calc.sgv(nblocks = nclusts, blk.sizes = obsperclust, vmat = SigHat)

  } else if('lme' %in% class(model)){

    X=stats::model.matrix(eval(model$call$fixed)[-2],
      data = model$data[,which( !(names(model$data)%in%c('zz','invwt')) )])
    n <- nrow(X)

    # Get grouping information from the model
    clust.id = names(summary(model)$groups)[1]
    obsperclust = as.numeric(table(model$data[,clust.id]))
    nclusts = length(obsperclust)

    SGV = calc.sgv(nblocks = nclusts, blk.sizes = obsperclust,
      vmat = bdiag(extract.lme.cov2(model, model$data, start.level=1)[['V']]))
  }

  # Use Helland's formula from Helland (1984)
  # with SGV instead of sigma^2 to form R^2_{SGV}

  SumSq = beta %*% var(X) %*% beta
  r2 = as.numeric(SumSq / (SumSq + SGV))

  if(adj==T){

    nprms = AICcmodavg::AICc(model, return.K = T)
    const = n - mean(obsperclust)
    r2adj = 1 - (1 - r2) * (const-1) / (const - nprms - 1)
    r2 = c('r2sgv' = r2, 'r2sgv.adj' = r2adj)

  }

  return(r2)

}






