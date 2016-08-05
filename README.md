# r2glmm

An R package for computation of model R squared and semi-partial 
R squared (with confidence limits) in linear and generalized linear mixed models

This package computes model and semi partial R squared and 
for the linear and generalized linear mixed model (LMM and GLMM). 
The R squared measure from Edwards et.al (2008) is extended
to the GLMM using penalized quasi-likelihood (PQL) estimation 
(see Jaeger et al. 2016). Three methods of computation are provided:

(1) The Kenward-Roger approach**,
(2) The method introduced by Nakagawa and Schielzeth (2013)****, and
(3) an approach using standardized generalized variance (SGV)
that can be used for both mean model and covariance model selection.

Confidence limits for semi partial R squared and model R squared are
computed for each of the methods listed.

Instructions for installation:

After installing the devtools package, run this code in the console:

devtools::install_github('bcjaeger/r2glmm')

**Due to some inconsistency between the pbkrtest package and the glmmPQL
function, the Kenward-Roger approach in the r2glmm package is limited to
the LMM.
****The r2glmm package only computes marginal R squared for the LMM and
does not generalize the statistic to the GLMM; however, confidence limits
may be computed for this marginal R squared in the LMM using this package

