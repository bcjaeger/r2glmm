

#' @export

r2beta.glmerMod <- function(model, partial=TRUE, method='sgv',
                            data = NULL){

  if(is.null(data)) data = model@frame

  fam = model@resp$family
  frm = as.list(model@call)[['formula']]
  y = all.vars(stats::update(frm, . ~ 1))

  # Remove missing values (glm does this automatically)
  data = data[stats::complete.cases(data[,y]), ]

  fit0 = stats::glm(lme4::nobars(frm), family = fam, data = data)

  eta <- fit0$linear.predictors

  data$zz <- eta + fit0$residuals
  wz=NULL; data$wz <- fit0$weights

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

  return(r2beta(model = fit, method = method, partial = partial))

}
