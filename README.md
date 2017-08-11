
<!-- README.md is generated from README.Rmd. Please edit that file -->
r2glmm
======

This package computes model and semi partial *R*<sup>2</sup> with confidence limits for the linear and generalized linear mixed model (LMM and GLMM). The *R*<sup>2</sup> measure from Edwards et.al (2008) is extended to the GLMM using penalized quasi-likelihood (PQL) estimation (see Jaeger et al. 2016).

-   Changes: Version 0.1.2

1.  Optimized computation of matrix inverses and cross-products in order to decrease computation time.

-   Changes: Version 0.1.1

1.  Updated the r2beta function with an optional data argument. Users who wish to use a function (i.e. log(x)) in their model formula should specify the original data frame used by the model when using the r2beta function.
2.  Included support for GLMMs fitted using the glmer function from the lme4 package.
3.  Added generic plot and print functions for R2 objects.

-   Why use this package?

The *R*<sup>2</sup> statistic is a well known tool that describes goodness-of-fit for a statistical model. In the linear model, *R*<sup>2</sup> is interpreted as the proportion of variance in the data explained by the fixed predictors and semi-partial *R*<sup>2</sup> provide standardized measures of effect size for subsets of fixed predictors. In the linear mixed model, numerous definitions of *R*<sup>2</sup> exist and interpretations vary by definition. The r2glmm package computes *R*<sup>2</sup> using three definitions:

1.  *R*<sub>*β*</sub><sup>2</sup>, a standardized measure of multivariate association between the fixed predictors and the observed outcome. This method was introduced by Edwards et al. (2008)
2.  *R*<sub>*Σ*</sub><sup>2</sup>, the proportion of generalized variance explained by the fixed predictors. This method was introduced by Jaeger et al. (2017)
3.  *R*<sub>(*m*)</sub><sup>2</sup>, the proportion of variance explained by the fixed predictors. This method was introduced by Nakagawa and Schielzeth (2013) and later modified by Johnson (2014).

Each interpretation can be used for model selection and is helpful for summarizing model goodness-of-fit. While the information criteria are useful tools for model selection, they do not quantify goodness-of-fit, making the *R*<sup>2</sup> statistic an excellent tool to accompany values of AIC and BIC. Additionally, in the context of mixed models, semi-partial *R*<sup>2</sup> and confidence limits are two useful and exclusive features of the r2glmm package.

-   Instructions for installation:

The most up-to-date version of the r2glmm package is available on Github. To download the package from Github, after installing and loading the devtools package, run the following code from the R console:

``` r
devtools::install_github('bcjaeger/r2glmm')
```

Alternatively, There is a version of the package available on CRAN. To download the package from CRAN, run the following code from the R console:

``` r
install.packages('r2glmm')
```

-   How to use this package

The main function in this package is called r2beta. The r2beta function summarizes a mixed model by computing the model *R*<sup>2</sup> statistic and semi-partial *R*<sup>2</sup> statistics for each fixed predictor in the model. The r2glmm package computes *R*<sup>2</sup> using three definitions. Below we list the methods, their interpretation, and an example of their application:

1.  *R*<sub>*β*</sub><sup>2</sup>, a standardized measure of multivariate association between the fixed predictors and the observed outcome. This statistic is primarily used to select fixed effects in the linear and generalized linear mixed model.

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

# Compute mean models with the r2beta statistic 
# using the Kenward-Roger approach.

m1 = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
m2 = lmer(distance ~ age + (1|Subject), data = Orthodont)

(r2.m1 = r2beta(m1, method = 'kr', partial = T))
#>    Effect   Rsq upper.CL lower.CL
#> 1   Model 0.671    0.771    0.563
#> 2     age 0.578    0.691    0.454
#> 4 age:Sex 0.074    0.212    0.004
#> 3     Sex 0.004    0.065    0.000
(r2.m2 = r2beta(m2, method = 'kr', partial = T))
#>   Effect   Rsq upper.CL lower.CL
#> 1  Model 0.589    0.698    0.468
#> 2    age 0.589    0.698    0.468
```

1.  *R*<sub>*Σ*</sub><sup>2</sup>, the proportion of generalized variance explained by the fixed predictors. This statistic is primarily used to select a covariance structure in the linear and generalized linear mixed model.

``` r

# m1 has a compound symmetric (CS) covariance structure.
m1 = lme(distance ~ age*Sex,  ~1|Subject, data = Orthodont)

# m2 is an order 1 autoregressive (AR1) model with
# gender-specific residual variance estimates.
m2 = lme(distance ~ age*Sex, data=Orthodont, 
         correlation = corAR1(form=~1|Subject),
         weights = varIdent(form=~1|Sex))

# Compare the models
(r2m1 = r2beta(m1,method='sgv'))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.559    0.669    0.447
#> 2           age 0.392    0.527    0.256
#> 4 age:SexFemale 0.038    0.144    0.000
#> 3     SexFemale 0.004    0.067    0.000
(r2m2 = r2beta(m2,method='sgv'))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.616    0.713    0.514
#> 2           age 0.454    0.580    0.323
#> 4 age:SexFemale 0.051    0.165    0.001
#> 3     SexFemale 0.006    0.075    0.000
```

1.  *R*<sub>(*m*)</sub><sup>2</sup>, the proportion of variance explained by the fixed predictors. This statistic is a simplified version of *R*<sub>*β*</sub><sup>2</sup> that can be used as a substitute for models fitted to very large datasets.

``` r

# Compute the R2 statistic using Nakagawa and Schielzeth's approach.
(r2nsj = r2beta(m1, method = 'nsj', partial = TRUE))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.410    0.540    0.290
#> 2           age 0.261    0.398    0.137
#> 4 age:SexFemale 0.021    0.105    0.000
#> 3     SexFemale 0.002    0.055    0.000

# Check the result with MuMIn's r.squaredGLMM
r2nsj_mum = MuMIn::r.squaredGLMM(m1)

all.equal(r2nsj[1,'Rsq'],as.numeric(r2nsj_mum[1]), tolerance = 1e-3)
#> [1] TRUE
```

-   *R*<sup>2</sup> for the Generalized Linear Mixed Model (GLMM)

The r2glmm package can compute *R*<sub>*β*</sub><sup>2</sup> for models fitted using the glmer function from the lme4 package. Note that this method is experimental in R and values of *R*<sub>*β*</sub><sup>2</sup> sometimes exceed 1. We recommend using the SAS macro available at <https://github.com/bcjaeger/R2FixedEffectsGLMM/blob/master/Glimmix_R2_V3.sas>. *R*<sub>*Σ*</sub><sup>2</sup> is more stable and can be computed for models fitted using either the glmer function or the glmmPQL function from the MASS package; however, minor differences in model estimation can lead to slight variation in the values of *R*<sub>*Σ*</sub><sup>2</sup>.

``` r

library(lattice)
library(MASS)
cbpp$period = as.numeric(cbpp$period)

# using glmer (based in lme4)
gm1 <- glmer(
  formula=cbind(incidence, size-incidence) ~ poly(period,2) + (1|herd),
  data = cbpp, family = binomial)

# using glmmPQL (based on nlme)
pql1 <- glmmPQL(
  cbind(incidence, size-incidence) ~ poly(period,2), 
  random = ~ 1|herd, family = binomial, data = cbpp
)
#> iteration 1
#> iteration 2
#> iteration 3
#> iteration 4

# Note minor differences in R^2_Sigma
r2beta(model = gm1, method = 'sgv', data = cbpp)
#>             Effect   Rsq upper.CL lower.CL
#> 1            Model 0.191    0.420    0.048
#> 2 poly(period, 2)1 0.180    0.397    0.029
#> 3 poly(period, 2)2 0.016    0.161    0.000
r2beta(model = pql1, method = 'sgv', data = cbpp)
#>             Effect   Rsq upper.CL lower.CL
#> 1            Model 0.210    0.438    0.059
#> 2 poly(period, 2)1 0.194    0.411    0.036
#> 3 poly(period, 2)2 0.025    0.182    0.000
```
