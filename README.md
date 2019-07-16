cvGEE: Cross-Validated Predictions from GEE
================

[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/cvGEE.svg?branch=master)](https://travis-ci.org/drizopoulos/cvGEE) [![CRAN status](http://www.r-pkg.org/badges/version/cvGEE)](https://cran.r-project.org/package=cvGEE) [![](https://cranlogs.r-pkg.org/badges/grand-total/cvGEE)](https://CRAN.R-project.org/package=cvGEE) [![Download counter](http://cranlogs.r-pkg.org/badges/cvGEE)](https://cran.r-project.org/package=cvGEE) 
[![Rdoc](http://www.rdocumentation.org/badges/version/cvGEE)](http://www.rdocumentation.org/packages/cvGEE)
[![rpackages.io rank](https://www.rpackages.io/badge/cvGEE.svg)](https://www.rpackages.io/package/cvGEE)

<img src="man/figures/logo.png" height="205" align="right"/>

Description
------------

<strong>cvGEE</strong> calculates cross-validated versions of the logarithmic, quadrative and spherical scoring rules for categorical data based on generalized estimating equations.

Basic Features
------------

- The package contains a single model-fitting function named `mixed_model()` with four 
required arguments, `fixed` a formula for the fixed effects, `random` a formula for the
random effects, `family` a family object specifying the type of response variable, and 
`data` a data frame containing the variables in the previously mentioned formulas.

- Methods for standard generics are provided, i.e., `coef()`, `fixef()`, `ranef()`, 
`vcov()`, `logLik()`, `summary()`, `anova()`, `confint()`, `fitted()`, `residuals()`, 
`predict()`, and `simulate()`.

- Negative binomial mixed models can be fitted using the `negative.binomial()` family 
object.

- Zero-inflated Poisson and negative binomial models using the `zi.poisson()` and 
`zi.negative.binomial()` family objects.

- Hurdle Poisson and negative binomial models using the `hurdle.poisson()` and 
`hurdle.negative.binomial()` family objects.

- Two-part/hurdle mixed models for semi-continuous normal data using the 
`hurdle.lognormal()` family objects.

- Continuation ratio mixed models for ordinal data using functions `cr_setup()` and `cr_marg_probs()`.

- Beta and hurdle Beta mixed effects models using `beta.fam()` and `hurdle.beta.fam()` 
family objects.

- Users may also specify their own log-density function for the repeated measurements 
response variable, and the internal algorithms will take care of the optimization.

- Calculates the marginalized coefficients using the idea of Hedeker et al. (2017) using 
function `marginal_coefs()`.

- Predictions with confidence interval for constructing effects plots are provided by 
function `effectPlotData()`.

Basic Use
------------

We compare a linear and a nonlinear GEE for the dichotomized version of serum bilirubin from the PBC dataset
```r
library("geepack")
library("splines")
library("lattice")

pbc2$serBilirD <- as.numeric(pbc2$serBilir > 1.2)

gm1 <- geeglm(serBilirD ~ year * drug, 
              family = binomial(), data = pbc2, id = id, 
              corstr = "exchangeable")

gm2 <- geeglm(serBilirD ~ ns(year, 3, Boundary.knots = c(0, 10)) * drug, 
              family = binomial(), data = pbc2, id = id, 
              corstr = "exchangeable")

plot_data <- cv_gee(gm1, return_data = TRUE)
plot_data$linear <- plot_data$.score
plot_data$non_linear <- unlist(cv_gee(gm2))

xyplot(linear + non_linear ~ year | .rule, data = plot_data, 
       type = "smooth", auto.key = TRUE, layout = c(3, 1),
       scales = list(y = list(relation = "free")),
       xlab = "Follow-up time (years)", ylab = "Scoring Rules")
```

Installation
------------

The development version of the package can be installed from GitHub using the **devtools**
package:
```r
devtools::install_github("drizopoulos/cvGEE")
```

and with vignettes
```r
devtools::install_github("drizopoulos/cvGEE", build_opts = NULL)
```
