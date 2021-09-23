# MEDdist: Multivariate Expectile-based Distribution (MED)

This package accompanies the article "Multivariate expectile-based distribution: properties, Bayesian inference and applications" by J. Arbel, S. Girard, H.D. Nguyen, and A. Usseglio-Carleve.

## Installation

The package can be installed from R via:

``` r
# install.packages("devtools")
devtools::install_github("AntoineUC/MEDdist")
```

## Usage

* Density, distribution function, quantile function and random generation for the MED distribution.
* Maximum likelihood estimation.
* Bayesian inference using Hamiltonian Monte Carlo (HMC) via RStan.

The following code samples from the MED and draws contour lines.

```{r}
library(MEDdist)
set.seed(1)

# sample from MED
mu <- c(0,0)
rho <- .8
N_iter <- 1000
MED_sample <- rmed(N_iter, mu = mu, rho = rho)

# scatter plot
plot(MED_sample, pch = 20, asp = 1, xlab = "", ylab = "", col = rgb(red = .1, green = .1, blue = .1, alpha = .3))
abline(h = 0, v = 0, col = "gray")
points(MED_sample)

# contour lines of MED
theta <- seq(0,2*pi, length = 1000)
level_vec <- seq(1/1000,.99, length = 10)
level_max <- max(level_vec)

for(level in level_vec){
  r <- sqrt(log(1/level)/(1+rho*cos(theta)))
  x <- r*cos(theta)
  y <- r*sin(theta)
  lines(x = x, y = y, lwd = 2, col = rgb(level/level_max, 0, (level_max-level)/level_max, alpha = .8))
}
```
