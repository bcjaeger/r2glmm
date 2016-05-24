

<<<<<<< HEAD
#' Compute r2beta with a specified C matrix
=======
#' Compute R2 with a specified C matrix
>>>>>>> 963695e99d0c3235e00058662fba7428a951e9d5
#' @param c Contrast matrix for fixed effects
#' @param x Fixed effects design matrix
#' @param SigHat estimated model covariance matrix
#' @param beta fixed effects estimates
#' @param NS.adj if True, the denominator degrees of freedom are
#' adjusted to replicate the R squared proposed by Nakagawa and Schielzeth.
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
cmp.R2 = function(c, x, SigHat, beta, NS.adj = F){

  # Compute relevant quantities for degrees of freedom
  rank.c = sum(apply(c, 2, sum))
  n = nrow(x)
  Xmt = matrix(x, nrow = n)

  # Compute the approximate Wald F Statistic
  if(NS.adj==T){
    denom = SigHat * as.matrix(c %*% solve( t(Xmt)%*%Xmt ) %*% t(c))
  } else {
    denom = as.matrix(c %*% solve(t(Xmt)%*%solve(SigHat)%*%Xmt ) %*% t(c))
  }

  wald = t(c %*% beta) %*% ginv(denom) %*% c%*%beta / rank.c

  # Compute the R2 statistic
  ddf = ifelse(NS.adj, n-1, n-1-rank.c)
  ndf = rank.c
  ss = ndf/ddf * wald
  R2 = ss/(1+ss)

  return(c(F = wald, v1 = ndf, v2 = ddf, ncp = wald*ndf, Rsq = R2))

}
