

#' @export


r2beta.lmerMod <- function(model, partial=TRUE, method){

    # Get model matrices
    X = lme4::getME(model, 'X')
    n <- nrow(X)

    # Get grouping information from the model
    clust.id = names(model@flist)[ length(model@flist) ]
    obsperclust = as.numeric(table(model@frame[ , clust.id ]))
    mobs = mean(obsperclust)
    nclusts = length(obsperclust)

    # The Kenward Roger Approach
    if (toupper(method) == 'KR') {

      ### Get random effects from model call

      random = paste('(',unlist(lme4::findbars(model@call)), ')',
                     sep = '',
                     collapse = '+')

      # Null model formula:
      # same covariance structure with fixed effects removed (except intercept)

      null.form <- stats::as.formula(paste('. ~ 1 +', random))

      # Compute Kenward Roger approximate F test using null model defined above

      mc <- pbkrtest::KRmodcomp(model,
               stats::update(model, null.form, data = model@frame))$stat

      # Compute the R2beta statistic for the full model
      # Store results in a dataframe

      R2 = data.frame(Effect = 'Model',
                      F = mc$Fstat, v1 = mc$ndf, v2 = mc$ddf,
                      ncp = mc$Fstat * mc$ndf,
                      Rsq = with(mc, (ndf*Fstat / ddf) / (1 + ndf*Fstat/ddf) ) )

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
        r2part = within(p,{
          Rsq = (p$num.Df*p$F/p$den.Df) / (1+p$num.Df*p$F/p$den.Df)
          ncp = p$F*p$num.Df
        })


        names(r2part)[names(r2part)=="num.Df"] <- "v1"
        names(r2part)[names(r2part)=="den.Df"] <- "v2"

        R2 = rbind(R2, r2part)

      }

    }
    # SGV and NSJ methods
    else if (toupper(method) == 'SGV' | toupper(method) == 'NSJ'){

      # Get fixed effects estimates
      beta = lme4::fixef(model)
      p <- length(beta)

      # Get random effects design matrix
      Z = as.matrix(lme4::getME(model, 'Z'))

      # Get variance component estimates
      s2e  = lme4::getME(model, 'sigma')^2
      lam  = lme4::getME(model, 'Lambda')
      lamt = lme4::getME(model, 'Lambdat')

      G = s2e * as.matrix(lam %*% lamt)

      # Compute estimated covariance matrix
      SigHat = Z%*%G%*%t(Z) + s2e*diag(nrow(Z))

      if(toupper(method)=='NSJ'){

        # NSJ approach takes the mean of the trace of the model covariance
        SigHat = mean(diag(SigHat))

      }

      if(toupper(method)=='SGV'){

        # SGV approach takes standardized determinant of the model covariance
        SigHat = calc_sgv(nblocks = nclusts, vmat = SigHat)

      }

      # the C matrix defines the Wald F Test for Fixed Effects
      C = list(); nms = c('Model', names(beta)[-1])
      C[['Model']] = cbind(rep(0, p-1),diag(p-1))

      if (partial == T){

        # Add the partial contrast matrices to C
        for(i in 2:(p)){
          C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
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

  return(R2)

}

