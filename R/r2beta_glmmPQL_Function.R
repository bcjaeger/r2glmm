
#' @export r2beta.glmmPQL

r2beta.glmmPQL <- function(model, partial = T, method='sgv'){

  if(!require("stats",quietly=T)) stop("package 'stats' is essential")
  if(!require("dplyr",quietly=T)) stop("package 'dplyr' is essential")
  if(!require("mgcv",quietly=T)) stop("package 'mgcv' is essential")
  if(!require("MASS",quietly=T)) stop("package 'MASS' is essential")

  # Get model matrices
  X=stats::model.matrix(eval(model$call$fixed)[-2],
      data = model$data[,which( !(names(model$data)%in%c('zz','invwt')) )])
  n <- nrow(X)

  # Get grouping information from the model
  clust.id = names(summary(model)$groups)[1]
  obsperclust = as.numeric(table(model$data[,clust.id]))
  mobs = mean(obsperclust)
  nclusts = length(obsperclust)

  if(toupper(method)!='SGV'){
    stop('Only the SGV method is compatible with glmmPQL objects')
  }

  # Get fixed effects
  beta = fixef(model)
  p <- length(beta)

  # Get covariance matrix from the model

  mlist = mgcv::extract.lme.cov2(model, model$data, start.level=1)[['V']]

  SigHat = calc.sgv(nblocks = nclusts, vmat = mlist)

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
  r2=lapply(C, FUN=cmp.R2, x=X, SigHat=SigHat, beta=beta, method=method,
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

