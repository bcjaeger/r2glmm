

#' Compute r2beta with a specified C matrix
#' Compute R2 with a specified C matrix
#' @param c Contrast matrix for fixed effects
#' @param x Fixed effects design matrix
#' @param SigHat estimated model covariance (matrix or scalar)
#' @param beta fixed effects estimates
#' @return A vector with the Wald statistic (ncp), approximate Wald F
#' statistic (F), numerator degrees of freedom (v1), denominator degrees
#' of freedom (v2), and the specified r squared value (Rsq)
#' @examples
#' library(nlme)
#' library(lme4)
#' lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#' X = model.matrix(lmemod, data = Orthodont)
#' SigHat = extract.lme.cov(lmemod, data = Orthodont)
#' beta = fixef(lmemod)
#' p = length(beta)
#'
#' C = cbind(rep(0, p-1),diag(p-1))
#' partial.c = make.partial.C(p-1,p,2)
#'
#' cmp.R2(c=C, x=X, SigHat=SigHat, beta=beta)
#' cmp.R2(c=partial.c, x=X, SigHat=SigHat, beta=beta)
#' @export cmp.R2

cmp.R2 = function(c, x, SigHat, beta, method, obsperclust, nclusts){

  scalar = !is.matrix(SigHat)

  # Compute relevant quantities for degrees of freedom
  rank.c = sum(apply(c, 2, sum))
  num.obs = nrow(x)
  Xmt = matrix(x, nrow = num.obs)

  # Compute the approximate Wald F Statistic
  if(scalar){
    denom = SigHat * as.matrix(c %*% Matrix::solve( crossprod(Xmt) ) %*% t(c))
  } else {
    denom = as.matrix(c %*% solve(t(Xmt)%*%Matrix::solve(SigHat)%*%Xmt) %*% t(c))
  }

  wald = t(c %*% beta) %*% MASS::ginv(denom) %*% c%*%beta / rank.c

  # Compute the R2 statistic

  if(toupper(method)=='NSJ') ddf = num.obs - 1

  if(toupper(method)=='SGV'){

    mobs = mean(obsperclust)
    m = nclusts - mobs - 1

    if(m <= 0){

      ddf = num.obs - 1

      print('Number of observations per cluster is greater than the number of clusters. Setting ddf to n-1.')
    }

    if( m > 0){
      ddf = mobs * (m+1) + (mobs-1)*(mobs-2) / 2
    }

  }

  ndf = rank.c
  ss = ndf/ddf * wald
  R2 = ss/(1+ss)

  return(c(F = wald, v1 = ndf, v2 = ddf, ncp = wald*ndf, Rsq = R2))

}


# Older version
# cmp.R2 = function(c, x, SigHat, beta){
#
#   scalar = !is.matrix(SigHat)
#
#   # Compute relevant quantities for degrees of freedom
#   rank.c = sum(apply(c, 2, sum))
#   num.obs = nrow(x)
#   Xmt = matrix(x, nrow = num.obs)
#
#   # Compute the approximate Wald F Statistic
#   if(scalar){
#     denom = SigHat * as.matrix(c %*% Matrix::solve( t(Xmt)%*%Xmt ) %*% t(c))
#   } else {
#     denom = as.matrix(c %*% solve(t(Xmt)%*%Matrix::solve(SigHat)%*%Xmt) %*% t(c))
#   }
#
#   wald = t(c %*% beta) %*% MASS::ginv(denom) %*% c%*%beta / rank.c
#
#   # Compute the R2 statistic
#   ddf = ifelse(scalar, num.obs-1, num.obs-1-rank.c)
#   ndf = rank.c
#   ss = ndf/ddf * wald
#   R2 = ss/(1+ss)
#
#   return(c(F = wald, v1 = ndf, v2 = ddf, ncp = wald*ndf, Rsq = R2))
#
# }
