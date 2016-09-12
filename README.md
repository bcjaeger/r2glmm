
<!-- README.md is generated from README.Rmd. Please edit that file -->
r2glmm
======

This package computes model and semi partial R squared with confidence limits for the linear and generalized linear mixed model (LMM and GLMM). The R squared measure from Edwards et.al (2008) is extended to the GLMM using penalized quasi-likelihood (PQL) estimation (see Jaeger et al. 2016).

-   Why use this package?

The \(R^2\) statistic is a well known tool that describes a statistical models goodness-of-fit. \(R^2\) may be interpreted as a measurement of the proportion of variance in the data explained by the fitted model. This quantity is often of interest for investigators in social and biological sciences. Since the measure is standardized, model \(R^2\) can be compared across studies for similar models. Thus, the \(R^2\) may be used for meta-analyses. Semi-partial \(R^2\) can be used to provide standardized measures of effect size for individual predictors, and allow investigators to select a set of predictors based on both statistical significance and relative importance.

Currently information criteria dominate the applied practice of selecting the most parsimonious mixed model. These criteria provide guidance, but cannot be used to measure goodness-of-fit. Further, they can only be used to compare models fitted to the same data. Lastly, the information criteria do not allow investigators to assess individual fixed predictors. Thus, it is beneficial to apply the information criteria in conjunction with \(R^2\) statistics when conducting statistical inference on mixed models.

-   Instructions for installation:

Currently, the r2glmm package is available at my github site. After installing and loading the devtools package, run this code from the R console:

devtools::install\_github('bcjaeger/r2glmm')

The files should then be downloaded and installed.

-   How to use this package

The main function in this package is called r2beta. A user may fit a mixed model for one of the supported model types and then apply the r2beta function using the specified model as input. Additionally, the investigator may specify whether semi-partial \(R^2\) are computed (they are by default) and what type of method to employ for computation. Three methods of computation are currently provided:

1.  an approach using standardized generalized variance (SGV) that can be used for covariance model selection. The SGV approach chooses

``` r

library(lme4)
#> Loading required package: Matrix
library(nlme)
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:lme4':
#> 
#>     lmList
library(r2glmm)

data(Orthodont)

# Linear mixed models
mermod = lmer(formula = distance ~ age*Sex + (1|Subject), data = Orthodont)

# Compute the R2 statistic using the SGV method
# This method and others are explained below
# Semi-partial R Squared are calculated by default
sgv_r2 = r2beta(mermod, method = 'sgv')

# Clean up the output a bit
nmrc = sapply(sgv_r2,is.numeric)
sgv_r2[,nmrc] = apply(sgv_r2[,nmrc], 2, round, 4)

sgv_r2
#>          Effect       F v1 v2      ncp    Rsq upper.CL lower.CL
#> 1         Model 40.1745  3 95 120.5234 0.5592   0.6691   0.4475
#> 2           age 61.1659  1 95  61.1659 0.3917   0.5267   0.2557
#> 4 age:SexFemale  3.7636  1 95   3.7636 0.0381   0.1435   0.0004
#> 3     SexFemale  0.3424  1 95   0.3424 0.0036   0.0665   0.0000
```

1.  The Kenward-Roger approach applies the small sample approximation to the F statistic using the pbkrtest package, and is recommended for selecting fixed effects. Due to some inconsistency between the pbkrtest package and the glmmPQL function, the Kenward-Roger approach in the r2glmm package is limited to the LMM.

``` r

# Compute the R2 statistic using the Kenward-Roger Approach.
(kr_r2 = r2beta(mermod, method = 'kr', partial = FALSE))
#>   Effect        F v1       v2      ncp       Rsq upper.CL  lower.CL
#> 1  Model 45.37291  3 66.64665 136.1187 0.6713115 0.770938 0.5631402
```

1.  The method introduced by Nakagawa and Schielzeth (2013) (see note B) computes the marginal R squared described by Nakagawa and Schielzeth, which was later modified by Johnson (2014). Additionally, this package computes confidence limits and semi-partial \(R^2\) for the marginal coefficient of determination. The r2glmm package only computes marginal R squared for the LMM and does not generalize the statistic to the GLMM.

``` r

# Compute the R2 statistic using the Kenward-Roger Approach.
nsj_r2 = r2beta(mermod, method = 'nsj', partial = TRUE)

# Clean up the output a bit
nmrc = sapply(nsj_r2,is.numeric)
nsj_r2[,nmrc] = apply(nsj_r2[,nmrc], 2, round, 4)

nsj_r2
#>          Effect       F v1  v2     ncp    Rsq upper.CL lower.CL
#> 1         Model 24.7691  3 107 74.3072 0.4098   0.5402   0.2901
#> 2           age 37.7111  1 107 37.7111 0.2606   0.3983   0.1366
#> 4 age:SexFemale  2.3204  1 107  2.3204 0.0212   0.1053   0.0001
#> 3     SexFemale  0.2111  1 107  0.2111 0.0020   0.0546   0.0000


# The MuMIn package computes the marginal R squared for fixed effects,
# but does not compute confidence limits or semi-partials

MuMIn::r.squaredGLMM(mermod)
#>       R2m       R2c 
#> 0.4098416 0.7827266
```
