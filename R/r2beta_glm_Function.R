
#' @export

r2beta.glm <- function(model, partial=TRUE, method='sgv', data = NULL){

  mod.pql = glmPQL(model)

  beta = stats::coef(mod.pql)
  p = length(beta)
  X = stats::model.matrix(mod.pql)
  SigHat = summary(mod.pql)$sigma^2

  # C matrix defines the Wald Test for Fixed Effects
  C = list(); nms = c('Model', names(beta)[-1])

  # Define the model Wald statistic for all fixed effects
  C[['Model']] = cbind(rep(0, p-1),diag(p-1))

  # For partial R2 statistics:
  if (partial == T){

    # add the partial contrast matrices to C
    for(i in 2:(p)) {
      C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
    }

  }

  # Compute the specified R2
  r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat, beta=beta, method='lm')

  # initialize a dataframe to hold results
  R2 = data.frame(Effect = names(r2))

  # place results in the dataframe
  for(i in names(r2[[1]])){
    R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
  }

  R2 = within(R2, {
            lower.CL = stats::qbeta(0.025, R2$v1/2, R2$v2/2, R2$ncp)
            upper.CL = stats::qbeta(0.975, R2$v1/2, R2$v2/2, R2$ncp)
            } )
  R2 = R2[order(-R2$Rsq),]

  class(R2) <- c('R2', 'data.frame')

  return(R2)

}
