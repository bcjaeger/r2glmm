

#'  R Squared Difference Test (R2DT). Test for a statistically significant difference in generalized explained variance between two candidate models.
#'
#' @param x An R2 object from the r2beta function.
#' @param y An R2 object from the r2beta function. If y is not specified, Ho: E[x] = mu is tested (mu is specified by the user).
#' @param mu Used to test Ho: E[x] = mu.
#' @param fancy if TRUE, the output values are rounded and changed to characters.
#' @param cor if TRUE, the R squared statistics are assumed to be positively correlated and a simulation based approach is used. If FALSE, the R squared are assumed independent and the difference of independent beta distributions is used. This only needs to be specified when two R squared measures are being considered.
#' @param onesided if TRUE, the alternative hypothesis is that one model explains a larger proportion of generalized variance. If false, the alternative is that the amount of generalized variance explained by the two candidate models is not equal.
#' @param nsims number of samples to draw when simulating correlated non-central beta random variables. This parameter is only relevant if cor=TRUE.
#' @param clim Desired confidence level for interval estimates regarding the difference in generalized explained variance.
#' @return A confidence interval for the difference in R Squared statistics and a p-value corresponding to the null hypothesis of no difference.
#' @examples
#' library(nlme)
#' library(lme4)
#' library(r2glmm)
#'
#' data(Orthodont)
#'
#' # Comparing two linear mixed models
#' m1 = lmer(distance ~ age*Sex+(1|Subject), Orthodont)
#' m2 = lmer(distance ~ age*Sex+(1+age|Subject), Orthodont)
#'
#' m1r2 = r2beta(model=m1,partial=FALSE)
#' m2r2 = r2beta(model=m2,partial=FALSE)
#'
#' # Accounting for correlation can make a substantial difference.
#'
#' r2dt(x=m1r2, y = m2r2, cor = TRUE)
#' r2dt(x=m1r2, y = m2r2, cor = FALSE)
#'
#'
#' @export

r2dt = function(x, y=NULL, cor = TRUE, fancy=FALSE,
                onesided=TRUE, clim=95,
                nsims=2000, mu=NULL){
  if(!('R2' %in% class(x) ))
    x = r2beta(x, partial = FALSE)
  if(!is.null(y) & !('R2'%in%class(y)))
    y = r2beta(y, partial = FALSE)


  alpha = ifelse(clim>1, 1-clim/100, 1-clim)
  limits = c(0+alpha/2, 1-alpha/2)

  pbetas <- function (z,a=1,b=1,c=1,d=1,ncp1=0,ncp2=0){
    n=length(z)
    tt=z
    for (i in 1:length(z)){
      f=function (xx) {
        dbetas(xx,a=a,b=b,c=c,d=d,ncp1=ncp1,ncp2=ncp2)
      }
      tt[i]=stats::integrate(f,lower=-1,upper=z[i])$value
    }
    return(tt)
  }

  dbetas <- function (z,a=1,b=1,c=1,d=1,ncp1=0,ncp2=0){
    n=length(z)
    tt=rep(0,n)
    for (i in 1:length(z))
    {f=function (xx) {
      stats::dbeta(xx,shape1=a,shape2=b,ncp=ncp1)*
        stats::dbeta(xx-z[i],shape1=c,shape2=d,ncp=ncp2)
    }
    if (z[i]>0&z[i]<=1){
      tt[i]=stats::integrate(f,lower=z[i],upper=1)$value
    }
    if (z[i]>=-1&z[i]<=0){
      tt[i]=stats::integrate(f,lower=0,upper=1+z[i])$value}
    }
    return(tt)
  }

  if(is.null(y)){

    # Testing one R squared value against a target

    if(is.null(mu)){
      stop('Specify mu in order to test E[R^2] = mu')
    }

    xx = x[1,]

    alpha  = xx$v1/2
    beta   = xx$v2/2
    lambda = xx$ncp

    xvals = seq(0,1, length.out = 10000)

    repeat{

      yvals = stats::dbeta(xvals, alpha, beta, lambda)

      E_r2 = mean(xvals*yvals)

      diff = E_r2-mu

      if(abs(diff)<0.0001) break

      # if mu < E[R^2]
      if (diff < 0){
        lambda=lambda*(1+abs(diff))
      } else {
        lambda=lambda*(1-diff)
      }

    }

    t = xx$Rsq

    # Ho: E[R^2] = mu vs. Ha: E[R^2] < mu
    if(t > mu){
      pval = ifelse(onesided, 1, 2) * (1-stats::pbeta(t, alpha, beta, lambda))
    } else {
      pval = ifelse(onesided, 1, 2) * (stats::pbeta(t, alpha, beta, lambda))
    }

    diff = t-mu

    res = list(d=diff,
               ci=c(xx[,'lower.CL'], xx[,'upper.CL']) - mu,
               p=pval)

  } else {

    xx = x[1,]; yy = y[1,]; diff = xx$Rsq - yy$Rsq

    if(cor){

      # It is assumed that there is a positive
      # correlation between R squared statistics
      # when they have identical mean models.

      rsq = min(xx$Rsq, yy$Rsq)
      ncp = min(xx$ncp, yy$ncp)

      beta_cor = 1 - (rsq / (1+rsq))^4

      samples = MASS::mvrnorm(
        n=nsims,
        mu = c(0,0),
        Sigma = matrix(c(1,beta_cor, beta_cor,1),
                       nrow=2, byrow=TRUE)
      )

      u = stats::pnorm(samples)

      df = data.frame(
        r2x = stats::qbeta(u[,1],
                           shape1 = yy$v1/2,
                           shape2 = yy$v2/2,
                           ncp    = ncp),
        r2y = stats::qbeta(u[,2],
                           shape1 = yy$v1/2,
                           shape2 = yy$v2/2,
                           ncp    = ncp))

      diffs = df$r2x-df$r2y
      diffs = diffs[!is.na(diffs)]

      if (diff < 0){
        pval = ifelse(onesided, 1, 2) * mean(diffs < diff)
      } else {
        pval = ifelse(onesided, 1, 2) * mean(diffs > diff)
      }

      #ci

      df = data.frame(
        r2x = stats::qbeta(u[,1],
                           shape1 = xx$v1/2,
                           shape2 = xx$v2/2,
                           ncp    = xx$ncp),
        r2y = stats::qbeta(u[,2],
                           shape1 = yy$v1/2,
                           shape2 = yy$v2/2,
                           ncp    = yy$ncp))

      diffs = df$r2x-df$r2y
      diffs = diffs[!is.na(diffs)]

      ci = as.numeric(stats::quantile(diffs,probs=limits))

    } else {

      # Null Distribution
      # The difference distribution centered on 0
      # Using the CDF, calculate P(r2diff > observed)
      # multiply by 2 if test is two sided


      # Find a dominant Rsq
      dom = eval(parse(
        text=names(which.max(c(xx = xx[,'Rsq'], yy = yy[,'Rsq'])))))

      pval <- ifelse(onesided, 1, 2) *
        (1 - pbetas(z = abs(diff),
                    a = dom$v1/2, b = dom$v2/2,
                    c = dom$v1/2, d = dom$v2/2,
                    ncp1 = dom$ncp, ncp2 = dom$ncp))

      r2_diff_ci = function(start){

        lower = upper = start

        # Upper confidence limit

        repeat{

          qq = pbetas(z = upper,
                      a = xx$v1/2, b = xx$v2/2,
                      c = yy$v1/2, d = yy$v2/2,
                      ncp1 = xx$ncp, ncp2 = yy$ncp)

          if (qq >= limits[2] | abs(qq-limits[2])<0.0001) break

          upper = upper+0.005

        }

        # Lower confidence limit

        repeat{

          qq = pbetas(z = lower,
                      a = xx$v1/2, b = xx$v2/2,
                      c = yy$v1/2, d = yy$v2/2,
                      ncp1 = xx$ncp, ncp2 = yy$ncp)

          if (qq <= limits[1] | abs(qq-limits[1])<0.0001) break

          lower = lower-0.005

        }

        return(c(lower, upper))

      }

      ci = r2_diff_ci(start = diff)

    }

    res = list(d=diff, ci=ci,p=pval)

  }

  if (fancy == T){

    make.ci = function(x, upper=NULL, lower=NULL, dig = 2){

      # make.ci gives nice confidence intervals

      fr = function(x, dig=2){

        # fr formats and rounds numbers.

        res=trimws(paste(format(round(x, dig), nsmall=dig)))
        return(res)
      }

      c1 = is.null(upper)
      c2 = is.null(lower)
      c3 = length(x)==3

      parens <- function(left, right = NULL){

        if(is.null(right)){

          paste0('(',fr(left),")")

        } else {

          paste0('(',fr(left),' ',fr(right),')')

        }

      }

      if(c1&c2&c3){

        ci = paste(fr(x[1]), parens(left=x[2], right = x[3]))

      } else {

        ci = paste(fr(x), parens(left=lower, right=upper))

      }

      return(ci)

    }

    res$p = ifelse(res$p < 0.001,
                   '<0.001',
                   sprintf("%.3f",round(pval,3)))
    res$ci = make.ci(diff, upper = res$ci[2], lower=res$ci[1])

  }


  return(res)

}




