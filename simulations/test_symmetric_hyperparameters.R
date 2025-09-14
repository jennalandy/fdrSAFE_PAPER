# remove.packages("fdrSAFE")
# remotes::install_github("jennalandy/fdrSAFE", lib = "~/apps/R_4.1.0/")

.libPaths("~/apps/R_4.1.0/")
source("test.R")

library(fdrSAFE)
library(MASS)
library(ggplot2)

pi0 = 0.8
sd0 = 1
lower = 1.33
upper = 4
prop_left = 1/2
n = 1000

null = function(n) {
  sample = rnorm(n, mean = 0, sd = sd0)
  return(sample)
}

alt = function(n) {
  n_left = floor(prop_left*n)
  sample = c(
    runif(n_left, -1*upper, -1*lower),
    runif(n - n_left, lower, upper)
  )
  return(sample)
}

true_fdr = function(test_statistics) {
  fdr = pi0 * dnorm(test_statistics, mean = 0, sd = sd0) / (
    pi0 * dnorm(test_statistics, mean = 0, sd = sd0) + 
    (1-pi0) * (
      prop_left * dunif(test_statistics, -1*upper, -1*lower) + 
      (1 - prop_left) * dunif(test_statistics, lower, upper)
    )
  )
  return(fdr)
}

empirical_Fdr = function(test_statistics, truth){
  # Fdr = proportion null among those assigned not-null
  abs_t = abs(test_statistics)
  Fdr = numeric(length(test_statistics))
  for (i in 1:length(test_statistics)) {
    t = test_statistics[i]
    if (sum(abs_t > abs(t)) == 0) {
      Fdr[i] = 0
    } else {
      Fdr[i] = sum((1-truth) * (abs_t > abs(t))) / sum(abs_t > abs(t))
    }
  }
  return(Fdr)
}

data_generating_fn = function() {
  dat = list(
    'zz' = c(null(round(n*pi0)), alt(round(n*(1-pi0)))),
    'truth' = c(rep(FALSE, round(n*pi0)), rep(TRUE, round(n*(1-pi0))))
  )
  dat$true_fdr = true_fdr(dat$zz)
  dat$true_Fdr = empirical_Fdr(dat$zz, dat$truth)
  return(dat)
}

to_pval_function = function(test_statistics) {
  sapply(test_statistics, function(t) {
    2*(1-pnorm(abs(t), 0, sd = sd0))
  })
}

run_inc_n_test(
  data_generating_fn = data_generating_fn, 
  N = 200, 
  path = "results/test_symmetric_hyperparameters",
  overwrite = TRUE,
  to_pval_function = function(test_statistics) {
    sapply(test_statistics, function(t) {
      2*(1-pnorm(abs(t), mean = 0, sd = sd0))
    })
  },
  type = 'asymmetric'
)
