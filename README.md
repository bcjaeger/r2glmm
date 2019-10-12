
<!-- README.md is generated from README.Rmd. Please edit that file -->

# r2glmm

This package computes model and semi partial \(R^2\) with confidence
limits for the linear and generalized linear mixed model (LMM and GLMM).
The \(R^2\) measure from Edwards et al. (2008) is extended to the GLMM
using penalized quasi-likelihood (PQL) estimation (see Jaeger et al.
(2016)).

  - Changes: Version 0.1.3

<!-- end list -->

1.  Updated semi-partial computation. Categorical variables and
    polynomial regression variables are now grouped.

<!-- end list -->

  - Changes: Version 0.1.2

<!-- end list -->

1.  Optimized computation of matrix inverses and cross-products in order
    to decrease computation time.

<!-- end list -->

  - Changes: Version 0.1.1

<!-- end list -->

1.  Updated the r2beta function with an optional data argument. Users
    who wish to use a function (i.e. log(x)) in their model formula
    should specify the original data frame used by the model when using
    the r2beta function.
2.  Included support for GLMMs fitted using the glmer function from the
    lme4 package.
3.  Added generic plot and print functions for R2 objects.

<!-- end list -->

  - Why use this package?

The \(R^2\) statistic is a well known tool that describes
goodness-of-fit for a statistical model. In the linear model, \(R^2\) is
interpreted as the proportion of variance in the data explained by the
fixed predictors and semi-partial \(R^2\) provide standardized measures
of effect size for subsets of fixed predictors. In the linear mixed
model, numerous definitions of \(R^2\) exist and interpretations vary by
definition. The r2glmm package computes \(R^2\) using three definitions:

1.  \(R_\beta^2\), a standardized measure of multivariate association
    between the fixed predictors and the observed outcome. This method
    was introduced by Edwards et al. (2008).
2.  \(R_\Sigma^2\), the proportion of generalized variance explained by
    the fixed predictors. This method was introduced by Jaeger et al.
    (2016)
3.  \(R_{(m)}^2\), the proportion of variance explained by the fixed
    predictors. This method was introduced by Nakagawa and Schielzeth
    (2013) and later modified by Johnson (2014).

Each interpretation can be used for model selection and is helpful for
summarizing model goodness-of-fit. While the information criteria are
useful tools for model selection, they do not quantify goodness-of-fit,
making the \(R^2\) statistic an excellent tool to accompany values of
AIC and BIC. Additionally, in the context of mixed models, semi-partial
\(R^2\) and confidence limits are two useful and exclusive features of
the r2glmm package.

  - Instructions for installation:

The most up-to-date version of the r2glmm package is available on
Github. To download the package from Github, after installing and
loading the devtools package, run the following code from the R console:

``` r
devtools::install_github('bcjaeger/r2glmm')
```

Alternatively, There is a version of the package available on CRAN. To
download the package from CRAN, run the following code from the R
console:

``` r
install.packages('r2glmm')
```

  - How to use this package

The main function in this package is called r2beta. The r2beta function
summarizes a mixed model by computing the model \(R^2\) statistic and
semi-partial \(R^2\) statistics for each fixed predictor in the model.
The r2glmm package computes \(R^2\) using three definitions. Below we
list the methods, their interpretation, and an example of their
application:

1.  \(R_\beta^2\), a standardized measure of multivariate association
    between the fixed predictors and the observed outcome. This
    statistic is primarily used to select fixed effects in the linear
    and generalized linear mixed model.

<!-- end list -->

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
library(splines)
data(Orthodont)

# Compute mean models with the r2beta statistic 
# using the Kenward-Roger approach.

mer1 = lmer(distance ~ bs(age)*Sex + (1|Subject), data = Orthodont)
mer2 = lmer(distance ~ age + (1|Subject), data = Orthodont)

(r2.mer1 = r2beta(mer1, method = 'kr', partial = T, data = Orthodont))
#>        Effect   Rsq upper.CL lower.CL
#> 1       Model 0.626    0.734    0.527
#> 2     bs(age) 0.586    0.702    0.465
#> 4 bs(age):Sex 0.086    0.256    0.021
#> 3         Sex 0.072    0.260    0.001
(r2.mer2 = r2beta(mer2, method = 'kr', partial = T, data = Orthodont))
#>   Effect   Rsq upper.CL lower.CL
#> 1  Model 0.589    0.698    0.468
#> 2    age 0.589    0.698    0.468
```

2.  \(R_\Sigma^2\), the proportion of generalized variance explained by
    the fixed predictors. This statistic is primarily used to select a
    covariance structure in the linear and generalized linear mixed
    model.

<!-- end list -->

``` r

# m1 has a compound symmetric (CS) covariance structure.

lme1 = lme(distance ~ age*Sex,  ~1|Subject, data = Orthodont)

# m2 is an order 1 autoregressive (AR1) model with
# gender-specific residual variance estimates.
lme2 = lme(distance ~ age*Sex, ~1|Subject, data=Orthodont, 
         correlation = corAR1(form=~1|Subject),
         weights = varIdent(form=~1|Sex))

# Compare the models
(r2m1 = r2beta(model=lme1,method='sgv',partial=FALSE))
#> Warning in model.matrix.default(~b$groups[[n.levels - i + 1]] - 1,
#> contrasts.arg = c("contr.treatment", : non-list contrasts argument ignored
#>   Effect   Rsq upper.CL lower.CL
#> 1  Model 0.559    0.669    0.447
(r2m2 = r2beta(model=lme2,method='sgv',partial=FALSE))
#> Warning in model.matrix.default(~b$groups[[n.levels - i + 1]] - 1,
#> contrasts.arg = c("contr.treatment", : non-list contrasts argument ignored
#>   Effect   Rsq upper.CL lower.CL
#> 1  Model 0.603    0.703    0.498
```

3.  \(R_{(m)}^2\), the proportion of variance explained by the fixed
    predictors. This statistic is a simplified version of \(R_\beta^2\)
    that can be used as a substitute for models fitted to very large
    datasets.

<!-- end list -->

``` r

# Compute the R2 statistic using Nakagawa and Schielzeth's approach.
(r2nsj = r2beta(mer1, method = 'nsj', partial = TRUE))
#>        Effect   Rsq upper.CL lower.CL
#> 1       Model 0.410    0.551    0.305
#> 2     bs(age) 0.263    0.409    0.149
#> 3         Sex 0.032    0.126    0.000
#> 4 bs(age):Sex 0.024    0.134    0.005

# Check the result with MuMIn's r.squaredGLMM
r2nsj_mum = MuMIn::r.squaredGLMM(mer1)
#> Registered S3 method overwritten by 'MuMIn':
#>   method         from
#>   predict.merMod lme4
#> Warning: 'r.squaredGLMM' now calculates a revised statistic. See the help
#> page.

all.equal(r2nsj[1,'Rsq'],as.numeric(r2nsj_mum[1]), tolerance = 1e-3)
#> [1] TRUE
```

  - \(R^2\) for the Generalized Linear Mixed Model (GLMM)

The r2glmm package can compute \(R_\beta^2\) for models fitted using the
glmer function from the lme4 package. Note that this method is
experimental in R and values of \(R_\beta^2\) sometimes exceed 1. We
recommend using the SAS macro available at
<https://github.com/bcjaeger/R2FixedEffectsGLMM/blob/master/Glimmix_R2_V3.sas>.
\(R_\Sigma^2\) is more stable and can be computed for models fitted
using either the glmer function or the glmmPQL function from the MASS
package; however, minor differences in model estimation can lead to
slight variation in the values of \(R_\Sigma^2\).

``` r

library(lattice)
library(MASS)
cbpp$period = as.numeric(cbpp$period)

# using glmer (based in lme4)
gm1 <- glmer(
  formula=cbind(incidence, size-incidence) ~ bs(period) + (1|herd),
  data = cbpp, family = binomial)

# using glmmPQL (based on nlme)
pql1 <- glmmPQL(
  cbind(incidence, size-incidence) ~ bs(period), 
  random = ~ 1|herd, family = binomial, data = cbpp
)
#> iteration 1
#> iteration 2
#> iteration 3
#> iteration 4

# Note minor differences in R^2_Sigma
r2beta(model = gm1, method = 'sgv', data = cbpp)
#>       Effect  Rsq upper.CL lower.CL
#> 1      Model 0.24    0.476    0.091
#> 2 bs(period) 0.24    0.476    0.091
r2beta(model = pql1, method = 'sgv', data = cbpp)
#> Warning in model.matrix.default(~b$groups[[n.levels - i + 1]] - 1,
#> contrasts.arg = c("contr.treatment", : non-list contrasts argument ignored
#>       Effect  Rsq upper.CL lower.CL
#> 1      Model 0.22    0.458    0.077
#> 2 bs(period) 0.22    0.458    0.077
```

# References

<div id="refs" class="references">

<div id="ref-edwards2008r2">

Edwards, Lloyd J, Keith E Muller, Russell D Wolfinger, Bahjat F Qaqish,
and Oliver Schabenberger. 2008. “An R2 Statistic for Fixed Effects in
the Linear Mixed Model.” *Statistics in Medicine* 27 (29): 6137–57.

</div>

<div id="ref-r2glmmJaeger">

Jaeger, Byron C., Lloyd J. Edwards, Kalyan Das, and Pranab K. Sen. 2016.
“An \(R^2\) statistic for fixed effects in the generalized linear mixed
model.” *Journal of Applied Statistics* 0 (0): 1–20.
<https://doi.org/10.1080/02664763.2016.1193725>.

</div>

<div id="ref-johnson2014extension">

Johnson, Paul CD. 2014. “Extension of Nakagawa & Schielzeth’s R2glmm to
Random Slopes Models.” *Methods in Ecology and Evolution* 5 (9): 944–46.

</div>

<div id="ref-nakagawa2013general">

Nakagawa, Shinichi, and Holger Schielzeth. 2013. “A General and Simple
Method for Obtaining R2 from Generalized Linear Mixed-Effects Models.”
*Methods in Ecology and Evolution* 4 (2): 133–42.

</div>

</div>
