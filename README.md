cvGEE: Cross-Validated Predictions from GEE
================

[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/cvGEE.svg?branch=master)](https://travis-ci.org/drizopoulos/cvGEE) [![CRAN status](http://www.r-pkg.org/badges/version/cvGEE)](https://cran.r-project.org/package=cvGEE) [![](https://cranlogs.r-pkg.org/badges/grand-total/cvGEE)](https://CRAN.R-project.org/package=cvGEE) [![Download counter](http://cranlogs.r-pkg.org/badges/cvGEE)](https://cran.r-project.org/package=cvGEE) 
[![Rdoc](http://www.rdocumentation.org/badges/version/cvGEE)](http://www.rdocumentation.org/packages/cvGEE)
[![rpackages.io rank](https://www.rpackages.io/badge/cvGEE.svg)](https://www.rpackages.io/package/cvGEE)

<img src="man/figures/logo.png" height="205" align="right"/>

Description
------------

<strong>cvGEE</strong> calculates cross-validated versions of the logarithmic, quadratic and spherical scoring rules based on generalized estimating equations.

The package presumes that the GEE has been solved using the `geeglm()` function of the [**geepack**](https://cran.r-project.org/package=geepack).

- For `family = gaussian()` only the quadratic rule is available calculated as the squared prediction error; lower values indicate a better predictive ability.

- For `family = binomial()` and `family = poisson()` the probabilities for each event/category are calculated using the binomial and Poisson probability mass functions, respectively. For these families all three scoring rules are available, with higher values in each rule indicating better predictive ability.


Basic Use
------------

We compare a linear and a nonlinear GEE for the dichotomized version of the serum bilirubin biomarker from the PBC dataset:
```r
library("geepack")
library("cvGEE")
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
