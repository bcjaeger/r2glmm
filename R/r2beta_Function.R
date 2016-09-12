
#------------------------------------------------------------------------------#
#' r2beta
#'
#' @description Computes coefficient of determination (R squared) from
#' edwards et al., 2008 and the generalized R squared from Jaeger et al., 2016.
#' Currently implemented for linear mixed models with
#' \code{\link{lmer}} and \code{\link{lme}} objects. For
#' generalized linear mixed models, only \code{\link{glmmPQL}} are supported.
#'
#' @param model a fitted mermod, lme, or glmmPQL model.
#'
#' @param partial  if TRUE, semi-partial R squared are calculated for each
#' fixed effect in the mixed model.
#'
#' @param method Specifies the method of computation for R squared beta:
#'            if \code{method} = 'sgv' then the standardized generalized variance
#'            approach is applied. This method is recommended for covariance model
#'            selection.
#'            if \code{method} = 'kr', then the Kenward Roger approach is applied.
#'            This option is only available for \code{\link{lme}} models.
#'            if \code{method} = 'nsj',then the Nakagawa and Schielzeth approach
#'            is applied. This option is available for
#'            \code{\link{lmer}} and \code{\link{lme}} objects.
#'            if \code{method} = 'lm', the classical R squared from the
#'            linear model is computed. This method should only be used
#'            on glm and lm object.
#' @return A dataframe containing the model F statistic, numerator
#'   and denominator degrees of freedom, non-centrality parameter,
#'   and R squared statistic with 95% confidence limits.
#'   If partial = TRUE, then the dataframe also contains partial R
#'   squared statistics for all fixed effects in the model.
#'
#' @references  Edwards, Lloyd J., et al. "An R2 statistic for fixed effects in
#'    the linear mixed model." Statistics in medicine 27.29 (2008): 6137-6157.
#'
#' Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for
#'    obtaining R2 from generalized linear mixed effects models." Methods in
#'    Ecology and Evolution 4.2 (2013): 133-142.
#'
#'Jaeger, Byron C., et al., "An R Squared Statistic for Fixed Effects in the
#'    Generalized Linear Mixed Model." Journal of Applied Statistics (2016).
#'
#' @examples
#' library(nlme)
#' library(lme4)
#' data(Orthodont)
#'
#' # Linear mixed models
#' mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
#' lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#'
#' # The Kenward-Roger approach
#' r2beta(mermod, method = 'kr')
#'
#' # Standardized Generalized Variance
#' r2beta(mermod, method = 'sgv')
#' r2beta(lmemod, method = 'sgv')
#'
#' # The marginal R squared by Nakagawa and Schielzeth (extended by Johnson)
#' r2beta(mermod, method = 'nsj')
#'
#' # linear and generalized linear models
#'
#' library(datasets)
#' dis = data.frame(discoveries)
#' dis$year = 1:nrow(dis)
#'
#' lmod = lm(discoveries ~ year + I(year^2), data = dis)
#' glmod = glm(discoveries ~ year + I(year^2), family = 'poisson', data = dis)
#'
#' # Using an inappropriate link function (normal) leads to
#' # a poor fit relative to the poisson link function.
#'
#' r2beta(lmod)
#' r2beta(glmod)
#'
#' # PQL models
#' # Currently only SGV method is supported
#' library(MASS)
#' PQL_bac = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#'                   family = binomial, data = bacteria,
#'                   verbose = FALSE)
#'
#' r2beta(PQL_bac, method='sgv')
#'
#' @rdname r2beta
#' @export r2beta
#------------------------------------------------------------------------------#

r2beta <- function(model, partial=TRUE, method){
  UseMethod('r2beta')
}



# Old Version -------------------------------------------------------------


# r2beta <- function(model, partial = T, method='kr'){
#
#   if(!require("stats",quietly=T)) stop("package 'stats' is essential")
#   if(!require("afex",quietly=T)) stop("package 'afex' is essential")
#   if(!require("pbkrtest",quietly=T)) stop("package 'pbkrtest' is essential")
#   if(!require("dplyr",quietly=T)) stop("package 'dplyr' is essential")
#   if(!require("mgcv",quietly=T)) stop("package 'mgcv' is essential")
#   if(!require("MASS",quietly=T)) stop("package 'MASS' is essential")
#
#   # If the model has no random effects (i.e. a linear model)
#
#   if (any(class(model) == 'lm')){
#     beta = coef(model)
#     sm = summary(model)
#     R2 = data.frame(Effect = 'Model', F = sm$fstatistic['value'],
#                     v1 = sm$fstatistic['numdf'],
#                     v2 = sm$fstatistic['dendf'],
#                     ncp = sm$fstatistic['value'] * sm$fstatistic['numdf'],
#                     Rsq = sm$r.squared)
#   }
#
#   # If the model is a mermod model (fit with lmer function)
#
#   else if (any(c('merModLmerTest', 'lmerMod') %in% class(model))){
#
#     # Get model matrices
#     X = getME(model, 'X')
#     n <- nrow(X)
#
#     # Get grouping information from the model
#     clust.id = names(model@flist)[ length(model@flist) ]
#     obsperclust = as.numeric(table(model@frame[,clust.id]))
#     mobs = mean(obsperclust)
#     nclusts = length(obsperclust)
#
#     # The Kenward Roger Approach
#     if (toupper(method) == 'KR') {
#
#       ### Get random effects from model call
#
#       random = paste('(',unlist(lme4::findbars(model@call)), ')',
#                      sep = '',
#                      collapse = '+')
#
#       # Null model formula:
#       # same covariance structure with fixed effects removed (except intercept)
#
#       null.form <- stats::as.formula(paste('. ~ 1 +', random))
#
#       # Compute Kenward Roger approximate F test using null model defined above
#
#       mc <- pbkrtest::KRmodcomp(model,
#                          update(model, null.form, data = model@frame))$stat
#
#       # Compute the R2beta statistic for the full model
#       # Store results in a dataframe r2
#
#       R2 = data.frame(Effect = 'Model',
#                       F = mc$Fstat, v1 = mc$ndf, v2 = mc$ddf,
#                       ncp = mc$Fstat * mc$ndf,
#                       Rsq = with(mc, (ndf*Fstat / ddf) / (1 + ndf*Fstat/ddf) ) )
#
#       # use mixed function to conduct KR F tests for each
#       # fixed effect in the full model
#
#       if (partial == T){
#
#         partials <- afex::mixed(model@call$formula,
#                                 data = model@frame,
#                                 progress = F)$anova_table
#
#         p = data.frame(partials)
#         p$Effect = rownames(partials)
#         p = p[, c('F', 'num.Df', 'den.Df', 'Effect')]
#
#         # Compute partial R2beta statistics for all fixed effects
#         r2part = dplyr::mutate(p,
#                   Rsq = (num.Df * F / den.Df) / (1 + num.Df*F/den.Df),
#                   ncp = F * num.Df
#                   ) %>% dplyr::rename(v1 = num.Df, v2 = den.Df)
#
#         R2 = rbind(R2, r2part)
#
#       }
#
#     }
#     # SGV and NSJ methods
#     else if (toupper(method) == 'SGV' | toupper(method) == 'NSJ'){
#
#       # Get fixed effects estimates
#       beta = fixef(model)
#       p <- length(beta)
#
#       # Get random effects design matrix
#       Z = as.matrix(getME(model, 'Z'))
#
#       # Get variance component estimates
#       s2e = getME(model, 'sigma')^2
#       G = s2e * as.matrix(getME(model, 'Lambda') %*% getME(model, 'Lambdat'))
#
#       # Compute estimated covariance matrix
#       SigHat = Z%*%G%*%t(Z) + s2e*diag(nrow(Z))
#
#       if(toupper(method)=='NSJ'){
#
#         # NSJ approach takes the mean of the trace of the model covariance
#         SigHat = mean(diag(SigHat))
#
#       }
#
#       if(toupper(method)=='SGV'){
#
#         # SGV approach takes standardized determinant of the model covariance
#         SigHat = calc.sgv(nblocks = nclusts, vmat = SigHat)
#
#       }
#
#       # the C matrix defines the Wald F Test for Fixed Effects
#       C = list(); nms = c('Model', names(beta)[-1])
#       C[['Model']] = cbind(rep(0, p-1),diag(p-1))
#
#       if (partial == T){
#
#         # Add the partial contrast matrices to C
#         for(i in 2:(p)){
#           C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
#         }
#
#       }
#
#       # Compute the specified R2
#       r2=lapply(C, FUN=cmp.R2, x=X, SigHat=SigHat, beta=beta, method=method,
#                 obsperclust=obsperclust, nclusts=nclusts)
#
#       # initialize a dataframe to hold results
#       R2 = data.frame(Effect = names(r2))
#
#       # place results in the dataframe
#       for(i in names(r2[[1]])){
#         R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
#       }
#     }
#   }
#
#   ### If the model is an lme object (fit with lme function)
#   else if ('lme' %in% class(model)){
#
#     # Get model matrices
#     X=stats::model.matrix(eval(model$call$fixed)[-2],
#            data = model$data[,which( !(names(model$data)%in%c('zz','invwt')) )])
#     n <- nrow(X)
#
#     # Get grouping information from the model
#     clust.id = names(summary(model)$groups)[1]
#     obsperclust = as.numeric(table(model$data[,clust.id]))
#     mobs = mean(obsperclust)
#     nclusts = length(obsperclust)
#
#     if('glmmPQL' %in% class(model) & toupper(method)=='NSJ'){
#     stop('r2beta only reduces exactly to the Nakagawa and Schielzeth R Squared in the Normal Linear Mixed Model
#          Using NS.adj = T on a glmmPQL object is not recommended')
#     }
#
#     # The Kenward Roger Approach can't be applied to these models
#     if (toupper(method) == 'KR'){
#       stop('The Kenward Roger approach is only compatible with merMod objects.')
#     }
#
#     if(toupper(method) == 'SGV' | toupper(method) == 'NSJ'){
#
#       # Get fixed effects
#       beta = fixef(model)
#       p <- length(beta)
#
#       # Get covariance matrix from the model
#
#       mlist = mgcv::extract.lme.cov2(model, model$data, start.level=1)[['V']]
#
#       if(toupper(method)=='NSJ'){
#
#         # NS approach takes the mean of the trace of SigHat
#         SigHat = mean(diag(as.matrix(bdiag(mlist))))
#
#       }
#
#       # SGV approach takes the standardized generalized variance of SigHat
#       if(toupper(method)=='SGV'){
#
#         SigHat = calc.sgv(nblocks = nclusts, vmat = mlist)
#
#       }
#
#       # C matrix defines the Wald Test for Fixed Effects
#       C = list(); nms = c('Model', names(beta)[-1])
#
#       # Define the model Wald statistic for all fixed effects
#       C[['Model']] = cbind(rep(0, p-1),diag(p-1))
#
#       # For partial R2 statistics:
#       if (partial == T){
#
#         # add the partial contrast matrices to C
#         for(i in 2:(p)) {
#           C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
#         }
#
#       }
#
#       # Compute the specified R2
#       r2=lapply(C, FUN=cmp.R2, x=X, SigHat=SigHat, beta=beta, method=method,
#                 obsperclust=obsperclust, nclusts=nclusts)
#
#       # initialize a dataframe to hold results
#       R2 = data.frame(Effect = names(r2))
#
#       # place results in the dataframe
#       for(i in names(r2[[1]])){
#         R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
#       }
#     }
#   }
#
#   R2 = mutate(R2,
#               lower.CL = qbeta(0.025, v1/2, v2/2, ncp),
#               upper.CL = qbeta(0.975, v1/2, v2/2, ncp)) %>%
#     dplyr::arrange(desc(Rsq))
#
#   return(R2)
#
# }
#
