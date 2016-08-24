
#' @export r2beta.lmerMod

r2beta.lmerMod <- function(model, partial = T, method='kr'){

  if(!require("stats",quietly=T)) stop("package 'stats' is essential")
  if(!require("afex",quietly=T)) stop("package 'afex' is essential")
  if(!require("pbkrtest",quietly=T)) stop("package 'pbkrtest' is essential")
  if(!require("dplyr",quietly=T)) stop("package 'dplyr' is essential")
  if(!require("MASS",quietly=T)) stop("package 'MASS' is essential")

    # Get model matrices
    X = getME(model, 'X')
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
                                update(model, null.form, data = model@frame))$stat

      # Compute the R2beta statistic for the full model
      # Store results in a dataframe r2

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
        r2part = dplyr::mutate(p, Rsq = (num.Df*F/den.Df) / (1+num.Df*F/den.Df),
                               ncp = F*num.Df) %>%
          dplyr::rename(v1 = num.Df, v2 = den.Df)

        R2 = rbind(R2, r2part)

      }

    }
    # SGV and NSJ methods
    else if (toupper(method) == 'SGV' | toupper(method) == 'NSJ'){

      # Get fixed effects estimates
      beta = fixef(model)
      p <- length(beta)

      # Get random effects design matrix
      Z = as.matrix(getME(model, 'Z'))

      # Get variance component estimates
      s2e = getME(model, 'sigma')^2
      G = s2e * as.matrix(getME(model, 'Lambda') %*% getME(model, 'Lambdat'))

      # Compute estimated covariance matrix
      SigHat = Z%*%G%*%t(Z) + s2e*diag(nrow(Z))

      if(toupper(method)=='NSJ'){

        # NSJ approach takes the mean of the trace of the model covariance
        SigHat = mean(diag(SigHat))

      }

      if(toupper(method)=='SGV'){

        # SGV approach takes standardized determinant of the model covariance
        SigHat = calc.sgv(nblocks = nclusts, vmat = SigHat)

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
      r2=lapply(C, FUN=cmp.R2, x=X, SigHat=SigHat, beta=beta, method=method,
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

  R2 = mutate(R2,
              lower.CL = qbeta(0.025, v1/2, v2/2, ncp),
              upper.CL = qbeta(0.975, v1/2, v2/2, ncp)) %>%
    dplyr::arrange(desc(Rsq))

  return(R2)

}

