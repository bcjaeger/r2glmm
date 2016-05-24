

#' Computes coefficient of determination (r squared) from edwards et al., 2008.
#'
#' @param model a fitted lmerMod object
#' @param data the data used for the fitted model.
#' @param partial  if TRUE, semi-partial R squared are calculated for each
#' fixed effect in the mixed model.
#' @param ddf Approximation method for the denominator degrees of freedom
#'            if ddf = 'res' then residual degrees of freedom are used
#'            if ddf = 'kr', then the Kenward Roger approach is applied.
#' @param NS.adj if TRUE, the denominator degrees of freedom are
#' adjusted to replicate the R squared proposed by Nakagawa and Schielzeth.
#' @return A dataframe containing the model F statistic, numerator
#'   and denominator degrees of freedom, and R^2 statistic. If
#'   partial = TRUE, then the dataframe also contains partial
#'   R^2 statistics for all fixed effects in the model
#' @references  Edwards, Lloyd J., et al. "An R2 statistic for fixed effects in
#'    the linear mixed model." Statistics in medicine 27.29 (2008): 6137-6157.
#'
#' Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for
#'    obtaining R2 from generalized linear mixed‐effects models." Methods in
#'    Ecology and Evolution 4.2 (2013): 133-142.
#' @examples
#' library(nlme)
#' library(lme4)
#' linmod = lm(distance~age*Sex, data = Orthodont)
#' mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
#' lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#' r2beta(linmod, data = Orthodont, partial = T)
#' r2beta(lmemod, data = Orthodont, partial = T, ddf = 'res')
#' r2beta(mermod, data = Orthodont, partial = T, ddf = 'res')
#' r2beta(mermod, data = Orthodont, partial = T, ddf = 'kr')
#'
#' # PQL models
#' # Currently only residual degrees of freedom supported
#' library(MASS)
#' PQL_mod = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#' family = binomial, data = bacteria, verbose = F)
#' r2beta(PQL_mod, data=bacteria, partial=T, ddf='res')
#' @export r2beta

r2beta <- function(model, data=NULL, partial = F, ddf = 'res', NS.adj = F){

  if(is.null(data)) stop('The data argument is necessary')
  if(NS.adj == T & ddf == 'kr'){
    ddf = 'res'; print('ddf=kr is not compatible with NS.adj is TRUE.')
  }
  if(!require("stats")) stop("package 'stats' is essential")
  if(!require("afex")) stop("package 'afex' is essential")
  if(!require("pbkrtest")) stop("package 'pbkrtest' is essential")
  if(!require("dplyr")) stop("package 'dplyr' is essential")
  if(!require("mgcv")) stop("package 'mgcv' is essential")
  if(!require("MASS")) stop("package 'MASS' is essential")


  ### If the model has no random effects (i.e. a linear model)

  if (any(class(model) == 'lm')){
    sm = summary(model)
    R2 = data.frame(Effect = 'Model', F = sm$fstatistic['value'],
                    v1 = sm$fstatistic['numdf'],
                    v2 = sm$fstatistic['dendf'],
                    ncp = sm$fstatistic['value'] * sm$fstatistic['numdf'],
                    Rsq = sm$r.squared)
  }

  ### If the model is a mermod model (fit with lmer function)

  else if (any(c('merModLmerTest', 'lmerMod') %in% class(model))){
    # The Kenward Roger Approach
    if (ddf %in% c('kr', 'KR')) {

      ### Get random effects from model call

      random = paste('(',unlist(lme4::findbars(model@call)), ')',
                     sep = '',
                     collapse = '+')

      # Null model formula:
      # same covariance structure with fixed effects removed except the intercept

      null.form <- stats::as.formula(paste('. ~ 1 +', random))

      ### Compute Kenward Roger approximate F test using null model defined above

      mc <- pbkrtest::KRmodcomp(model, update(model, null.form))$stat

      # Compute the R2beta statistic for the full model
      # Store results in a dataframe r2

      R2 = data.frame(Effect = 'Model',
                      F = mc$Fstat,
                      v1 = mc$ndf,
                      v2 = mc$ddf,
                      ncp = mc$Fstat * mc$ndf,
                      Rsq = with(mc, (ndf*Fstat / ddf) / (1 + ndf*Fstat / ddf) ) )

      # use mixed function to conduct KR F tests for each
      # fixed effect in the full model

      if (partial == T){

        partials <- afex::mixed(model@call$formula,
                                data = model@frame,
                                progress = F)$anova_table

        p = data.frame(partials)
        p$Effect = rownames(partials)
        p = p[, c('F', 'num.Df', 'den.Df', 'Effect')]

        # Compute partial R2beta statistics for all fixed effects

        r2part = dplyr::mutate(p, Rsq = (num.Df * F / den.Df) / (1 + num.Df*F/den.Df),
                               ncp = F * num.Df) %>% dplyr::rename(v1 = num.Df, v2 = den.Df)

        R2 = rbind(R2, r2part) %>% dplyr::arrange(desc(Rsq))

      }

    }
    # Residual Degrees of Freedom
    else if (ddf %in% c('res', 'Res')){

      # Get fixed effects
      beta = fixef(model)
      p <- length(beta)

      # Get model matrices
      X = stats::model.matrix(model, data = data)
      n <- nrow(X)
      Z = as.matrix(getME(model, 'Z'))

      # Covariance Matrix from the model
      s2e = getME(model, 'sigma')^2
      G = as.matrix(s2e * getME(model, 'Lambda') %*% getME(model, 'Lambdat'))
      SigHat = Z%*%G%*%t(Z) + s2e*diag(nrow(Z))

      # NS approach takes the mean of the trace of the model covariance
      if(NS.adj==T) SigHat = mean(diag(SigHat))

      # C matrix defines the Wald Test for Fixed Effects
      C = list(); nms = c('Model', names(beta)[-1])

      # Define the model Wald statistic for all fixed effects
      C[['Model']] = cbind(rep(0, p-1),diag(p-1))

      # For partial R2 statistics:
      if (partial == T){

        # add the partial contrast matrices to C
        for(i in 2:(p)){
          C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
        }

      }

      # Compute the specified R2
      r2=lapply(C, FUN=cmp.R2, x=X, SigHat, beta, NS.adj=NS.adj)

      # initialize a dataframe to hold results
      R2 = data.frame(Effect = names(r2))

      # place results in the dataframe
      for(i in names(r2[[1]])){
        R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
      }

      # Sort from highest R2 value to lowest
      # This is only meaningful if partial R2 are computed
      R2 = dplyr::arrange(R2, desc(Rsq))
    }
  }

  ### If the model is an lme object (fit with lme function)
  else if ('lme' %in% class(model)){

    if('glmmPQL' %in% class(model) & NS.adj==T){
    stop('r2beta only reduces exactly to the Nakagawa and Schielzeth R Squared in the Normal Linear Mixed Model
         Using NS.adj = T on a glmmPQL object is not recommended')
    }

    # The Kenward Roger Approach can't be applied to these models
    if (ddf %in% c('kr', 'KR')){
      stop('The Kenward Roger approach is only compatible with merMod objects.')
    }

    # Residual Degrees of Freedom
    if(ddf %in% c('res', 'Res')){

      # Get fixed effects
      beta = fixef(model)
      p <- length(beta)

      # Get model matrices
      X = stats::model.matrix(model, data = data)
      n <- nrow(X)

      # Get grouping information from the model
      cluster.name = names(summary(model)$groups)[1]
      obspersub = as.numeric(table(model$data[,cluster.name]))
      nsubs = length(obspersub)

      # Get covariance matrix from the model
      SigHat = bdiag(extract.lme.cov2(model, model$data, start.level=1)[['V']])

      # NS approach takes the mean of the trace of the model covariance
      if(NS.adj==T) SigHat = mean(diag(as.matrix(SigHat)))

      # C matrix defines the Wald Test for Fixed Effects
      C = list(); nms = c('Model', names(beta)[-1])

      # Define the model Wald statistic for all fixed effects
      C[['Model']] = cbind(rep(0, p-1),diag(p-1))

      # For partial R2 statistics:
      if (partial == T){

        # add the partial contrast matrices to C
        for(i in 2:(p)) C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)

      }

      # Compute the specified R2
      r2=lapply(C, FUN=cmp.R2, x=X, SigHat, beta, NS.adj=NS.adj)

      # initialize a dataframe to hold results
      R2 = data.frame(Effect = names(r2))

      # place results in the dataframe
      for(i in names(r2[[1]])){
        R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
      }

      # Sort from highest R2 value to lowest
      # This is only meaningful if partial R2 are computed
      R2 = dplyr::arrange(R2, desc(Rsq))

    }
  }



  return(R2)

}


