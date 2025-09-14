.libPaths("~/apps/R_4.1.0/")

source('helpers.R')
load("data/1_platinum_data.RData")
library(fdrSAFE)

library(locfdr)
library(fdrtool)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggdist)

color_list = list(
  "fdrSAFE" = "#E69F00",
  "locfdr" = "#D55E00",
  "fdrtool" = "#009E73",
  "qvalue" = "#0072B2"
)

results_dir = "results"

set.seed(321)
N = length(platinum_data$statistics)
B = 1000
all_pi0_metrics = matrix(nrow = 0, ncol = 6)
all_fdr_metrics = matrix(nrow = 0, ncol = 6)
all_calib_dat = matrix(nrow = 0, ncol = 13)
for (b in 1:B) {
  b_idx = sample(1:N, N, replace = TRUE)
  t_statistics = platinum_data$statistics[b_idx]
  hypothesis_labels = platinum_data$fold_change$DE[b_idx]
  Fdr = platinum_data$Fdr[b_idx]

  fdrSAFE_res <- NULL
  tryCatch({
    locfdr_res <- locfdr(t_statistics, pct0 = 0.001)
    fdrtool_res <- fdrtool(t_statistics, plot = 0)
    qvalue_res <- qvalue(p_from_t(t_statistics, df = 4))

    fdrSAFE_res <- fdrSAFE(
      t_statistics,
      df = 4,
      parallel = TRUE,
      type = 'asymmetric'
    )
  }, error = function(e) {
    print(e)
  })
  if (is.null(fdrSAFE_res)) {
    next
  }

  calib_dat <- get_calib_dat(
    fdrs = list(
      "fdrSAFE" = fdrSAFE_res$fdr,
      "locfdr" = locfdr_res$fdr,
      "fdrtool" = fdrtool_res$lfdr,
      "qvalue" = qvalue_res$lfdr
    ),
    truth = hypothesis_labels
  ) %>%
    ungroup() %>%
    dplyr::select(method, fdr_group, prop_null) %>%
    pivot_wider(
      names_from = fdr_group,
      values_from = prop_null
    )
  all_calib_dat <- rbind(all_calib_dat, calib_dat)

  pi0_metrics = list(
    'b' = b,
    'true' = mean(1 - hypothesis_labels),
    'fdrSAFE' = fdrSAFE_res$pi0,
    'locfdr' = unlist(locfdr_res$fp0['mlest','p0']),
    'fdrtool' = unname(fdrtool_res$param[,'eta0']),
    'qvalue' = qvalue_res$pi0
  )
  all_pi0_metrics = rbind(all_pi0_metrics, pi0_metrics)

  fdr_metrics = do.call(rbind, list(
    method_metrics(
      'fdrSAFE',
      estimated_fdr = fdrSAFE_res$fdr, 
      test_statistics = t_statistics,
      hypothesis_labels = hypothesis_labels,
      true_Fdr = Fdr
    ),
    method_metrics(
      'locfdr',
      estimated_fdr = locfdr_res$fdr, 
      test_statistics = t_statistics,
      hypothesis_labels = hypothesis_labels,
      true_Fdr = Fdr
    ),
    method_metrics(
      'fdrtool',
      estimated_fdr = fdrtool_res$lfdr, 
      test_statistics = t_statistics,
      hypothesis_labels = hypothesis_labels,
      true_Fdr = Fdr
    ),
    method_metrics(
      'qvalue',
      estimated_fdr = qvalue_res$lfdr, 
      test_statistics = t_statistics,
      hypothesis_labels = hypothesis_labels,
      true_Fdr = Fdr
    )
  ))
  fdr_metrics = cbind(fdr_metrics, rep(b, 4))
  all_fdr_metrics <- rbind(
    all_fdr_metrics, fdr_metrics
  )

  write.csv(all_calib_dat, file.path(results_dir, "bootstrap_calib_metrics.csv"))
  write.csv(all_fdr_metrics, file.path(results_dir, "bootstrap_fdr_metrics.csv"))
  write.csv(all_pi0_metrics, file.path(results_dir, "bootstrap_pi0_metrics.csv"))
}