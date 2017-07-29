
#' @export

r2beta.gls <- function(model, partial=TRUE, method='sgv',
                       wtdbin = FALSE, data = NULL){


  # The Kenward Roger Approach can't be applied to these models
  if (toupper(method) == 'KR'){
    stop('The Kenward Roger approach is only compatible with lmerMod objects.')
  }


  if(is.null(data)) data = eval(model$call$data)

  # Get model matrices
  X=stats::model.matrix(
    stats::as.formula(model$call[[2]]), data = data)

  # Get number of observations
  n <- nrow(X)

  # Get grouping information from the model

  if(!is.null(model$groups)){
    clust.id = levels(model$groups)
    obsperclust = as.numeric(table(model$groups))
    mobs = mean(obsperclust)
    nclusts = length(obsperclust)
  } else {
    stop('Only functional for gls models with group specification')
  }

    # Get fixed effects
    beta = model$coefficients
    p <- length(beta)

    # Get covariance matrix from the model

    mlist <- nlme::corMatrix(model$modelStruct$corStruct)

    if (!is.null(model$modelStruct$varStruct)) {
      vw <- 1 / nlme::varWeights(model$modelStruct$varStruct)
    } else {
      vw <- 1
    }

    vars <- (model$sigma * vw)^2

    if(length(vars)>1){

      stop = cumsum(obsperclust)
      start = c(1, stop[1: (length(stop)-1)] + 1)

      for(i in seq(nclusts)){

        sdMat = sqrt( diag(vars[start[i]:stop[i]]) )

        mlist[[i]] = sdMat %*% mlist[[i]] %*% t(sdMat)
      }

    } else {
      mlist = lapply(mlist, function(x) x*vars)
    }

    if(toupper(method)=='NSJ'){

      # NS approach takes the mean of the trace of SigHat
      SigHat = mean(diag(as.matrix(Matrix::bdiag(mlist))))

    }

    # SGV approach takes the standardized generalized variance of SigHat
    if(toupper(method)=='SGV'){

      SigHat = calc_sgv(nblocks = nclusts, vmat = mlist)

    }

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
    r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat, beta=beta,
              method=method, obsperclust=obsperclust, nclusts=nclusts)

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
