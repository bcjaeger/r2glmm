

#' Compute R2 with a specified C matrix
#' @param c Contrast matrix for fixed effects
#' @param x Fixed effects design matrix
#' @param SigHat estimated model covariance (matrix or scalar)
#' @param beta fixed effects estimates
#' @param method the method for computing r2beta
#' @param obsperclust number of observations per cluster (i.e. subject)
#' @param nclusts number of clusters (i.e. subjects)
#' @return A vector with the Wald statistic (ncp), approximate Wald F
#' statistic (F), numerator degrees of freedom (v1), denominator degrees
#' of freedom (v2), and the specified r squared value (Rsq)
#' @examples
#' library(nlme)
#' library(lme4)
#' library(mgcv)
#' lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#' X = model.matrix(lmemod, data = Orthodont)
#' SigHat = extract.lme.cov(lmemod, data = Orthodont)
#' beta = fixef(lmemod)
#' p = length(beta)
#' obsperclust = as.numeric(table(lmemod$data[,'Subject']))
#' nclusts = length(obsperclust)
#' C = cbind(rep(0, p-1),diag(p-1))
#' partial.c = make.partial.C(p-1,p,2)
#'
#' cmp_R2(c=C, x=X, SigHat=SigHat, beta=beta, obsperclust = obsperclust,
#' nclusts = nclusts, method = 'sgv')
#' cmp_R2(c=partial.c, x=X, SigHat=SigHat, beta=beta, obsperclust = obsperclust,
#' nclusts = nclusts, method = 'sgv')
#' @export cmp_R2

cmp_R2 = function(c, x, SigHat, beta, method,
                  obsperclust=NULL, nclusts=NULL){

  scalar = !is.matrix(SigHat)

  # Compute relevant quantities for degrees of freedom
  rank.c = sum(apply(c, 2, sum))
  num.obs = nrow(x)
  Xmt = Matrix::Matrix(x, sparse = TRUE)

  # Compute the R2 statistic

  if(toupper(method)=='LM') ddf = num.obs - length(beta)

  if(toupper(method)=='NSJ') ddf = num.obs - 1

  if(toupper(method)=='SGV'){

    mobs = mean(obsperclust)
    m = nclusts - mobs - 1
    if(m <= 0){
      ddf = num.obs - 1
    } else {
      ddf = mobs * (m+1) + (mobs-1)*(mobs-2) / 2
    }

  }

  # Compute the approximate Wald F Statistic
  if(scalar){

    XtXinv <- Matrix::solve(Matrix::crossprod(Xmt))
    denom = SigHat * (c%*%XtXinv%*%t(c))

  } else {

    XtSX_inv <- Matrix::solve(Matrix::t(Xmt)%*%Matrix::solve(SigHat)%*%Xmt)
    denom = (c%*% XtSX_inv %*%t(c))

  }

  wald = suppressWarnings(try(as.numeric(
    t(c %*% beta) %*% Matrix::solve(denom) %*% c%*%beta / rank.c),
    silent = TRUE
  ))

  if(class(wald)=='try-error'){
    wald = try(as.numeric(
      t(c %*% beta) %*% MASS::ginv(as.matrix(denom)) %*% c%*%beta / rank.c)
    )
  }

  if(class(wald)=='try-error')
    stop('Could not invert adjusted covariance of regression estimates')

  ndf = rank.c
  ss = ndf/ddf * wald
  R2 = ss/(1+ss)

  c(F = wald,v1 = ndf, v2 = ddf, ncp = wald*ndf, Rsq = R2)

}


