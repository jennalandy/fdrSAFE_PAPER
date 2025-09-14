
library(PRROC)
library(pROC)
library(dplyr)

nrow_null0 <- function(df) {
  if (is.null(df)) {
    n = 0
  } else {
    n = nrow(df)
  }
  return(n)
}

source("https://raw.githubusercontent.com/jennalandy/gridsemblefdr/master/R/run_fdrtool_row.R")
source("https://raw.githubusercontent.com/jennalandy/gridsemblefdr/master/R/run_locfdr_row.R")
source("https://raw.githubusercontent.com/jennalandy/gridsemblefdr/master/R/run_qvalue_row.R")
source("https://raw.githubusercontent.com/jennalandy/gridsemblefdr/master/R/run_row.R")
source("https://raw.githubusercontent.com/jennalandy/gridsemblefdr/master/R/utilities.R")

get_oracle_grid <- function(dat, to_pval_function, grids, topq) {
  method_list = c(
    rep('locfdr', nrow_null0(grids$locfdr)),
    rep('fdrtool', nrow_null0(grids$fdrtool)),
    rep('qvalue', nrow_null0(grids$qvalue))
  )
  row_list = unlist(lapply(
    unique(method_list),
    function(m) {
      seq_len(sum(method_list == m))
    }
  ))
  
  oracle_grid <- do.call(rbind, lapply(1:length(row_list), function(i) {
    row_res <- run_row(
      test_statistics = dat$zz,
      to_pval_function = to_pval_function,
      grids = grids,
      method = method_list[i],
      row = row_list[i]
    )
    if (!is.null(row_res)) {
      # if not null, record pi0 estimate and metrics
      row_cutoff <- quantile(row_res$fdr, 1 - min(1, row_res$pi0))
      this_metrics <- calc_metrics(
        test_statistics = dat$zz,
        fdr = row_res$fdr, 
        Fdr = row_res$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = row_cutoff
      )
      this_metrics$pi0 <- row_res$pi0
    } else {
      # o.w. record placeholder (NA) metrics and pi0 estimate
      this_metrics <- rep(NA, 17)
      names(this_metrics) <- c(
        "pr", "roc", "brier", "Fdrerror", "equal_class_Fdrerror",
        "pr_topq", "roc_topq", "brier_toq", 
        "Fdrerror_topq", "pred_pos", "cutoff",
        "accuracy", "precision", "recall", "specificity",
        "f1", "pi0"
      )
    }
    
    # record method and grid row
    this_metrics$method <- method_list[i]
    this_metrics$row <- row_list[i]
    
    return(data.frame(this_metrics))
  }))
  
  return(oracle_grid)
}


calc_metrics <- function(test_statistics, fdr, Fdr, truth, true_Fdr, topq, true_fdr = NULL, cutoff = 0.2) {
  pred = rep(1, length(test_statistics))
  pred[fdr > cutoff] = 0
  
  TP = sum(pred == 1 & truth == 1)
  FP = sum(pred == 1 & truth == 0)
  TN = sum(pred == 0 & truth == 0)
  FN = sum(pred == 0 & truth == 1)
  
  accuracy = (TP + TN)/(TP + TN + FP + FN)
  precision = TP/(TP+FP)
  recall = TP/(TP+FN) # fnr = 1-recall, same as sensitivity
  specificity = TN/(TN + FP) # fpr = 1-specificity
  f1 = (2 * precision * recall) / (precision+recall)
  
  null_Fdrerror <- get_MSE(estimate = Fdr[truth == 0], true = true_Fdr[truth == 0])
  alt_Fdrerror <- get_MSE(estimate = Fdr[truth == 1], true = true_Fdr[truth == 1])
  weights_Fdr <- 1 - 4*true_Fdr*(1-true_Fdr)

  metrics <- list(
    'pr' = get_prauc(fdr, truth),
    'roc' = get_roc(fdr, truth),
    'brier' = get_brier(fdr, truth),
    'Fdrerror' = get_MSE(Fdr, true_Fdr),
    'equal_class_Fdrerror' = mean(c(null_Fdrerror, alt_Fdrerror)),
    'confidence_weighted_Fdrerror' = get_MSE(Fdr, true_Fdr, weights_Fdr),
    'pr_topq' = get_prauc(fdr[topq], truth[topq]),
    'roc_topq' = get_roc(fdr[topq], truth[topq]),
    'brier_toq' = get_brier(fdr[topq], truth[topq]),
    'Fdrerror_topq' = get_MSE(Fdr[topq], true_Fdr[topq]),
    'pred_pos' = sum(pred == 1),
    'cutoff' = cutoff,
    'accuracy' = accuracy,
    'precision' = precision,
    'recall' = recall,
    'specificity' = specificity,
    'f1' = f1
  )
  
  if (!is.null(true_fdr)) {
    metrics$fdrerror <- get_MSE(fdr, true_fdr)
    metrics$fdrerror_topq <- get_MSE(fdr[topq], true_fdr[topq])

    null_fdrerror <- get_MSE(estimate = fdr[truth == 0], true = true_fdr[truth == 0])
    alt_fdrerror <- get_MSE(estimate = fdr[truth == 1], true = true_fdr[truth == 1])
    metrics$equal_class_fdrerror <- mean(c(null_fdrerror, alt_fdrerror))

    weights_fdr <- 1 - 4*true_fdr*(1-true_fdr)
    metrics$confidence_weighted_fdrerror <- get_MSE(fdr, true_fdr, weights_fdr)
  }
  
  return(metrics)
}


get_prauc <- function(fdr, truth) {
  tryCatch({
    pr <- PRROC::pr.curve(
      scores.class0 = fdr[!as.logical(truth)],
      scores.class1 = fdr[as.logical(truth)]
    )
    return(as.numeric(pr$auc.integral))
  }, error = function(e) {
    return(NA)
  })
}

get_roc <- function(fdr, truth) {
  
  if (sum(truth) == 0 | sum(!truth) == 0) {
    return(0)
  }
  
  # direction = ">" accounts for inverse
  # relationship between fdr and y
  tryCatch({
    r <- pROC::roc(
      truth ~ fdr,
      direction = ">",
      levels = levels(as.factor(truth))
    )
    return(as.numeric(r$auc))
  }, error = function(e) {
    return(NA)
  })
}

get_brier <- function(fdr, truth) {
  prob_1 = 1-fdr
  mean((prob_1 - as.numeric(truth))**2)
}

get_MSE <- function(estimate, true, weights = rep(1, length(true))) {
  sum(weights*((true - estimate)**2))/sum(weights)
}