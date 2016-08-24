
#' @export r2beta.lm

r2beta.lm <- function(model, partial = T, method = NULL){

    beta = coef(model)
    p = length(beta)
    X = model.matrix(model)
    SigHat = summary(model)$sigma^2

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
    r2=lapply(C, FUN=cmp.R2, x=X, SigHat=SigHat, beta=beta, method='lm',
              obsperclust=obsperclust, nclusts=nclusts)

    # initialize a dataframe to hold results
    R2 = data.frame(Effect = names(r2))

    # place results in the dataframe
    for(i in names(r2[[1]])){
      R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
    }

    R2 = mutate(R2,
                lower.CL = qbeta(0.025, v1/2, v2/2, ncp),
                upper.CL = qbeta(0.975, v1/2, v2/2, ncp)) %>%
      dplyr::arrange(desc(Rsq))

    return(R2)

}
