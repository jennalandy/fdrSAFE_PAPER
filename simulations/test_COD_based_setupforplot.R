# remove.packages("gridsemblefdr")
# install.packages("../gridsemblefdr_asym_working_model", repos = NULL, type = "source", lib = "~/apps/R_4.1.0/")

.libPaths("~/apps/R_4.1.0/")
library(gridsemblefdr)
library(curatedOvarianData)
source("test.R")
library(MASS)
library(ggplot2)

pi0 = 0.8
prop_alt_over = 0.8
meandiff = 2
sddiff = 0.5
y = c(rep(0, 10), rep(1, 10))

data(TCGA_eset)
X_dat <- Biobase::exprs(TCGA_eset)
X_dat_variable <- X_dat[apply(X_dat, 1, var) >= quantile(apply(X_dat, 1, var), 0.9),]
n = nrow(X_dat_variable)
Sigma = cov(t(X_dat_variable))
mu0 = rowMeans(X_dat_variable)

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
  n_over = round(n*(1-pi0)*prop_alt_over)
  n_under = round(n*(1-pi0)) - n_over
  diff_idx = sample(1:n, n_over + n_under)
  diff = 1:n %in% diff_idx
  over_idx = sample(diff_idx, n_over)
  over = 1:n %in% over_idx
  under_idx = setdiff(diff_idx, over_idx)
  under = 1:n %in% under_idx

  mu1 = mu0
  mu1[over] = mu1[over] + rnorm(sum(over), mean = meandiff, sd = sddiff)
  mu1[under] = mu1[under] - rnorm(sum(under), mean = meandiff, sd = sddiff)

  X = matrix(nrow = n, ncol = length(y))
  X[,y == 0] = t(mvrnorm(sum(y == 0), mu = mu0, Sigma = Sigma))
  X[,y == 1] = t(mvrnorm(sum(y == 1), mu = mu1, Sigma = Sigma))

  test_statistics = sapply(
    1:nrow(X), 
    function(i) {
      test = t.test(X[i,y == 1], X[i,y == 0])
      return(test$statistic)
    }
  )

  dat = list(
    'zz' = test_statistics,
    'truth' = diff
  )
  dat$true_Fdr = empirical_Fdr(dat$zz, dat$truth)
  return(dat)
}
