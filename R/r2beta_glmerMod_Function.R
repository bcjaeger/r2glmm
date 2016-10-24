

#' @export

r2beta.glmerMod <- function(model,  partial=TRUE, method='sgv'){

  fam = model@resp$family
  data = model@frame
  frm = as.list(model@call)[['formula']]
  wz=NULL

  fit0 = stats::glm(lme4::nobars(frm), family = fam, data = data)

  eta <- fit0$linear.predictors

  data$zz <- eta + fit0$residuals
  data$wz <- fit0$weights

  niter = 40

  for (i in seq_len(niter)) {

    if(i == niter) stop('Model failed to converge')

    fit <- lme4::lmer(formula = stats::update(frm, zz ~ .),
                      data = data, weights = wz)

    etaold <- eta
    eta    <- stats::fitted(fit)

    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) break

    mu <- fam$linkinv(eta)
    mu.eta.val <- fam$mu.eta(eta)

    data$zz <- eta + (fit0$y - mu)/mu.eta.val
    data$wz <- mu.eta.val^2/fam$variance(mu)

  }

  return(r2beta(fit, method = method, partial = partial))

}
