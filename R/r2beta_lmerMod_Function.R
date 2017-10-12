

#' @export


r2beta.lmerMod <- function(model, partial=TRUE, method='sgv',
                           data = NULL){

  if(is.null(data)) data = model@frame
  if(is.null(data) & partial)
    stop('Please specify the dataframe used to fit the model.')

  # Get model matrices
  X = lme4::getME(model, 'X')
  n <- nrow(X)

  # Get grouping information from the model
  clust.id = names(model@flist)[ length(model@flist) ]

  if(!clust.id%in%names(data)){

    ids = strsplit(clust.id, ':')[[1]]
    data[[clust.id]] = interaction(data[[ids[1]]], data[[ids[2]]])
    data = data = droplevels(data)

  }

  obsperclust = as.numeric(table(data[ , clust.id ]))
  mobs = mean(obsperclust)
  nclusts = length(obsperclust)

  # The Kenward Roger Approach
  if (toupper(method) == 'KR') {

    # Calculate F statistic using the approach of Kenward & Roger

    # random effects (re)
    re = paste(sapply(lme4::findbars(stats::formula(model)), function(x){
      x = as.character(x)
      paste0('(',x[2], x[1], x[3],')')
    }), collapse = '+')

    null_model = stats::as.formula(paste0('. ~ 1 + ', re))

    # model comparison (mc)
    mc = pbkrtest::KRmodcomp(model, stats::update(model, null_model))$stats

    ss = with(mc, ndf*Fstat/ddf)

    R2 = data.frame(Effect='Model', 'F' = mc$Fstat,
                    v1 = mc$ndf, v2 = mc$ddf,
                    pval = mc$p.value,
                    ncp = mc$ndf*mc$Fstat,
                    Rsq = ss / (1+ss))

    # For partial R2 statistics:
    if(partial == T){

      suppressMessages(
        expr = p <- afex::mixed(stats::formula(model),
                                data=data,
                                progress = FALSE)$anova_table
      )

      names(p)<-c('v1','v2','F','pval')
      ss = p$v1*p$F/p$v2

      R2 = rbind(R2,
                 data.frame(Effect = rownames(p),
                            'F' = p$F,
                            v1 = p$v1, v2 = p$v2,
                            pval = p$pval,
                            ncp = p$v1*p$F,
                            Rsq = ss / (1+ss)))

    }

  } else if (toupper(method) == 'SGV' | toupper(method) == 'NSJ'){

    beta = lme4::fixef(model)
    p <- length(beta)

    if(p==1) stop('Model must have at least one fixed effect')

    # Get random effects design matrix
    Z = lme4::getME(model, 'Z')

    # Get variance component estimates
    s2e  = lme4::getME(model, 'sigma')^2
    lam  = lme4::getME(model, 'Lambda')
    lamt = lme4::getME(model, 'Lambdat')

    G = s2e * (lam %*% lamt)

    # Compute estimated covariance matrix
    SigHat = Z %*% ( G%*%Matrix::t(Z) )

    # Add the residual component
    Matrix::diag(SigHat) = Matrix::diag(SigHat) + s2e/model@resp$weights

    if(toupper(method)=='NSJ'){

      # NSJ approach takes the mean of the trace of the model covariance
      SigHat = mean(Matrix::diag(SigHat))

    }

    if(toupper(method)=='SGV'){

      # SGV approach takes standardized
      # determinant of the model covariance
      SigHat=calc_sgv(nblocks=nclusts,
                      blksizes=obsperclust,
                      vmat=SigHat)

    }

    # C matrix defines the Wald Test for Fixed Effects
    C = list(); nms = c('Model', names(beta)[-1])

    # Define the model Wald statistic for all fixed effects

    C[['Model']] = cbind(rep(0, p-1),diag(p-1))

    # For partial R2 statistics:
    if (partial == T & p>1){

      asgn = attr(X, 'assign')
      nmrs = 1:length(asgn)
      assign = split(nmrs, asgn)
      nTerms = length(assign)
      labs = attr(stats::terms(model), 'term.labels')
      nms = c('Model', labs)
      names(assign) = c('(Intercept)',labs)

      # add the partial contrast matrices to C
      for(i in 2:(nTerms)) {
        C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = assign[[i]])
      }

    }

    # Compute the specified R2
    r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat, beta=beta, method=method,
              obsperclust=obsperclust, nclusts=nclusts)

    # initialize a dataframe to hold results
    R2 = data.frame(Effect = names(r2))

    # place results in the dataframe
    for(i in names(r2[[1]])){
      R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
    }

  }

  # Calculate confidence limits with the non-central beta quantile function.
  # Arrange the resulting dataframe from highest to lowest R^2 values.

  R2 = within(R2, {
    lower.CL = stats::qbeta(0.025, R2$v1/2, R2$v2/2, R2$ncp)
    upper.CL = stats::qbeta(0.975, R2$v1/2, R2$v2/2, R2$ncp)
  } )
  R2 = R2[order(-R2$Rsq),]

  class(R2) <- c('R2', 'data.frame')

  return(R2)

}
