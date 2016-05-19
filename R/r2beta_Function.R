
#' Computes R^2 statistic from edwards et al., 2008.
#'
#' @param model a fitted lmerMod object
#' @param partial  if TRUE, semi-partial R squared are calculated for each fixed effect in the mixed model.
#' @param ddf Approximation method for the denominator degrees of freedom
#'            if ddf = 'res' then residual degrees of freedom are used
#'            if ddf = 'kr', then the Kenward Roger approach is applied.
#' @return A dataframe containing the model F statistic, numerator
#'   and denominator degrees of freedom, and R^2 statistic. If
#'   partial = TRUE, then the dataframe also contains partial
#'   R^2 statistics for all fixed effects in the model
#' @examples
#' library(nlme)
#' library(lme4)
#' m = lmer(distance ~ age*Sex+ (1|Subject), data = Orthodont)
#' r2beta(m, partial = T, ddf = 'kr')
#' @export r2beta

r2beta <- function(model, partial = F, ddf = 'kr'){

  ### If the model has no random effects (i.e. a linear model)

  if (class(model) == 'lm'){
    sm = summary(model)
    R2 = data.frame(Effect = 'Model', F = sm$fstatistic['value'],
                    v1 = sm$fstatistic['numdf'],
                    v2 = sm$fstatistic['dendf'],
                    ncp = sm$fstatistic['value'] * sm$fstatistic['numdf'],
                    Rsq = sm$r.squared)
  }
  else if (ddf %in% c('kr', 'KR')) {

  ### Get random effects from model call

  random = paste('(',unlist(lme4::findbars(model@call)), ')',
                 sep = '',
                 collapse = '+')

  # Null model formula:
  # same covariance structure with fixed effects removed except the intercept

  null.form <- as.formula(paste('. ~ 1 +', random))

  ### Compute Kenward Roger approximate F test using null model defined above

  mc <- pbkrtest::KRmodcomp(model, update(model, null.form))$stat

  # Compute the R2beta statistic for the full model
  # Store results in a dataframe r2

  R2 = data.frame(Effect = 'Model',
                  F = mc$Fstat,
                  v1 = mc$ndf,
                  v2 = mc$ddf,
                  ncp = mc$Fstat * mc$ndf,
                  Rsq = with(mc, (ndf*Fstat / ddf) / (1 + ndf*Fstat / ddf) ) )

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

    r2part = dplyr::mutate(p, Rsq = (num.Df * F / den.Df) / (1 + num.Df*F/den.Df),
                    ncp = F * num.Df) %>% dplyr::rename(v1 = num.Df, v2 = den.Df)

    R2 = rbind(R2, r2part) %>% dplyr::arrange(desc(Rsq))

    }

  }
  else if (ddf %in% c('res', 'Res')){

    # Get fixed effects
    beta = fixef(model)
    p <- length(beta)

    # Get model matrices
    X = stats::model.matrix(model)
    n <- nrow(X)
    Z = getME(model, 'Z')

    # C matrix defines the Wald Test for Fixed Effects
    C = list(); nms = c('Model', names(beta)[-1])

    # Define the usual Wald statistic for all fixed effects
    C[['Model']] = cbind(rep(0, p-1),diag(p-1))

    G.list = list()

    for(part in names(model@flist)){

      # get the counts of units within blocks
      N = as.numeric(table(model@flist[part]))

      # take the covariance matrix for the current variance component
      g = matrix(VarCorr(model)[[part]], nrow = dim(VarCorr(model)[[part]])[1])

      # stack it, and then store it in the G.list
      G.list[[part]] = kronecker(diag(length(N)), g)

    }

    # Stacked (diagonal) random effects covariance matrix
    G = as.matrix(bdiag(G.list))

    # Covariance Matrix from the model
    SigHat = Z%*%G%*%t(Z) + sigma(model)^2*diag(nrow(Z))

    # For partial R2 statistics...
    if (partial == T){
      # define a function to generate partial contrast matrices
      make.partial.C = function(rows, cols, index){
        x=matrix(0, nrow = rows, ncol = cols)
        x[index-1, index] = 1
        return(x)
      }
      # add the partial contrast matrices to C
      for(i in 2:(p)){
        C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = i)
      }
    }

    # Define a function to compute R2 with a specified C matrix
    cmp.R2 = function(c){

      # Compute relevant quantities for degrees of freedom
      rank.c = sum(apply(c, 2, sum))
      X = matrix(X, nrow = n)

      # Compute the approximate Wald F Statistic
      denom = as.matrix(c %*% solve(t(X)%*%solve(SigHat)%*%X ) %*% t(c))
      wald = t(c %*% beta) %*% ginv(denom) %*% c%*%beta / rank.c

      # Compute the R2 statistic
      ddf = n-1-rank.c; ndf = rank.c; ss = ndf/ddf * wald
      R2 = ss/(1+ss)

      return(c(F = wald, v1 = ndf, v2 = ddf, ncp = wald*ndf, Rsq = R2))

    }

    # Compute the specified R2
    r2=lapply(C, cmp.R2)

    # initialize a dataframe to hold results
    R2 = data.frame(Effect = names(r2))

    # place results in the dataframe
    for(i in names(r2[[1]])){
      R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
    }

    # Sort from highest R2 value to lowest
    # This is only meaningful if partial R2 are computed
    R2 = dplyr::arrange(R2, desc(Rsq))
  }

  return(R2)

}


