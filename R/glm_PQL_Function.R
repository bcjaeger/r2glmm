
#' Compute PQL estimates for fixed effects from a generalized linear model.
#' @param glm.mod a generalized linear model fitted with the glm function.
#' @param niter maximum number of iterations allowed in the PQL algorithm.
#' @param data The data used by the fitted model. This argument is required
#'        for models with special expressions in their formula, such as
#'        offset, log, cbind(sucesses, trials), etc.
#' @return A glmPQL object (i.e. a linear model using pseudo outcomes).
#' @examples
#'
#' # Load the datasets package for example code
#' library(datasets)
#' library(dplyr)
#'
#' # We'll model the number of world changing discoveries per year for the
#' # last 100 years as a poisson outcome. First, we set up the data
#'
#' dat = data.frame(discoveries) %>% mutate(year = 1:length(discoveries))
#'
#' # Fit the GLM with a poisson link function
#' mod <- glm(discoveries~year+I(year^2), family = 'poisson', data = dat)
#'
#' # Find PQL estimates using the original GLM
#' mod.pql = glmPQL(mod)
#'
#' # Note that the PQL model yields a higher R Squared statistic
#' # than the fit of a strictly linear model. This is attributed
#' # to correctly modelling the distribution of outcomes and then
#' # linearizing the model to measure goodness of fit, rather than
#' # simply fitting a linear model
#'
#' summary(mod.pql)
#' summary(linfit <- lm(discoveries~year+I(year^2), data = dat))
#'
#' r2beta(mod.pql)
#' r2beta(linfit)
#'
#' @export glmPQL

glmPQL <- function(glm.mod, niter = 20, data = NULL){

  if (!c('GLM')%in%toupper(class(glm.mod))){
    stop('only models with class GLM are supported')
  }

  fit0 <- glm.mod

  if(is.null(data)){
    pql.dat = data.frame(fit0$data)
  } else {
    pql.dat = data
  }

  mod.call <- fit0$call
  mod.form=mod.call[['formula']]
  w <- fit0$prior.weights
  off <- mod.call[['offset']]
  if(is.null(off)) off = 0
  eta <- fit0$linear.predictors
  zz <- eta + fit0$residuals
  wz <- fit0$weights
  fam <- fit0$family

  mod.names <- names(mod.call)[-1L]

  keep <- is.element(mod.names,
        c("fixed", "data", "subset", "na.action", "control"))

  for (i in mod.names[!keep]) mod.call[[i]] <- NULL

  mod.form[[2L]] <- quote(zz)
  mod.call[["formula"]] <- mod.form
  mod.call[[1L]] <- quote(lm)
  pql.dat$zz <- zz
  pql.dat$invwt <- 1/wz
  mod.call$data <- pql.dat

  for (i in seq_len(niter)) {
    fit <- eval(mod.call)
    eta.old <- eta
    eta <- stats::fitted(fit) + off
    if (sum((eta - eta.old)^2) < 1e-06 * sum(eta^2))
      break
    mu <- fam$linkinv(eta)
    mu.eta.val <- fam$mu.eta(eta)
    pql.dat$zz <- eta + (fit0$y - mu)/mu.eta.val - off
    wz <- w * mu.eta.val^2/fam$variance(mu)
    pql.dat$invwt <- 1/wz
    mod.call$data <- pql.dat
  }

  attributes(fit$logLik) <- NULL
  fit$call <- fit0$call
  fit$family <- fit0$family
  fit$logLik <- as.numeric(NA)
  fit$data = pql.dat
  oldClass(fit) <- c("glmPQL", oldClass(fit))
  return(fit)
}

