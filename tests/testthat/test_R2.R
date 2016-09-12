

library(r2glmm)
context("Computing r2beta")

library(nlme)
library(lme4)
data(Orthodont)

# Linear mixed models
mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)

lme.r2 = r2beta(model = lmemod, method = 'sgv', partial = FALSE)
mer.r2 = r2beta(model = mermod, method = 'sgv', partial = FALSE)

r2.diff = lme.r2$Rsq - mer.r2$Rsq

test_that("r2beta works", {
  expect_true(r2.diff < 0.02)
})

