





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

  return(SGV)

}


#' Computes R squared statistic for a generalized linear mixed model using
#' penalized quasi likelihood (PQL) estimation and standardized generalized
#' variance (SGV)
#
#' @param model a fitted mixed model (lmerMod or lme object)
#' @param adj if TRUE, an adjusted R squared for model selection is provided.
#' @param data the dataframe used to fit the mixed model.
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
#' r2.sgv(PQL_mod, data=bacteria)
#' @export r2.sgv
#'

r2.sgv <- function(model, data, adj = F){

  # Get fixed effects
  beta = fixef(model)
  p <- length(beta)

  X = stats::model.matrix(model, data = data)
  n <- nrow(X)

  if(any(c('merModLmerTest', 'lmerMod') %in% class(model))){

    Z = as.matrix(getME(model, 'Z'))
    G.list = N = num.clusters = list()

    for(part in names(model@flist)){

      # get the counts of units within blocks
      N[[part]] = as.numeric(table(eval(model@call$data)[, part]))
      num.clusters[[part]] = length(N[[part]])

      # take the covariance matrix for the current variance component
      g = matrix(VarCorr(model)[[part]], nrow = dim(VarCorr(model)[[part]])[1])

      # stack it, and then store it in the G.list
      G.list[[part]] = kronecker(diag(num.clusters[[part]]), g)

    }

    # Stacked (diagonal) random effects covariance matrix
    G = as.matrix(bdiag(G.list))

    # Covariance Matrix from the model
    v.y = Z%*%G%*%t(Z) + sigma(model)^2*diag(nrow(Z))

    # Number of clusters in the covariance matrix of responses
    blks = unlist(num.clusters)
    obspersub = N[[names(which.min(blks))]]
    nsubs = min(blks)

    SGV = calc.sgv(nblocks = nsubs, blk.sizes = obspersub, vmat = v.y)

  } else if('lme' %in% class(model)){

    # Get grouping information from the model
    cluster.name = names(summary(model)$groups)[1]
    obspersub = as.numeric(table(model$data[,cluster.name]))
    nsubs = length(obspersub)

    require(mgcv)

    SGV = calc.sgv(nblocks = nsubs, blk.sizes = obspersub,
      vmat = bdiag(extract.lme.cov2(model, model$data, start.level=1)[['V']]))
  }

  # Use Helland's formula from Helland (1984)
  # with SGV instead of sigma^2 to form R^2_{SGV}

  SumSq = beta %*% var(X) %*% beta
  r2 = as.numeric(SumSq / (SumSq + SGV))

  if(adj==T){

    # load AICc package to count the number of parameters in the model
    require(AICcmodavg)

    nprms = AICc(model, return.K = T)
    const = n / (mean(obspersub))
    r2adj = 1 - (1 - r2) * (const-1) / (const - nprms - 1)
    r2 = c('r2sgv' = r2, 'r2sgv.adj' = r2adj)

  }

  return(r2)

}


