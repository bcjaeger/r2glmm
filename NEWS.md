## version 0.1.3

* Added a `NEWS.md` file to track changes to the package.

* Updated semi-partial computation. Categorical variables and polynomial regression variables are now grouped. 

## version 0.1.2

* Optimized computation of matrix inverses and cross-products in order to decrease computation time.

## version 0.1.1

* Updated the r2beta function with an optional data argument. Users who wish to use a function (i.e. log(x)) in their model formula should specify the original data frame used by the model when using the r2beta function.

* Included support for GLMMs fitted using the glmer function from the lme4 package.

* Added generic plot and print functions for R2 objects.

