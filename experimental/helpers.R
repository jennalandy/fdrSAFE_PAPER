library(dplyr)

get_subset <- function(pi0, n = 1000, data = platinum_data) {
  # indices of null (0) and not-null (1)
  idx = 1:length(data$fold_change$DE)
  idx_0 = idx[data$fold_change$DE == 0]
  idx_1 = idx[data$fold_change$DE == 1]

  n0 = round(n*pi0)
  n1 = n - n0

  # sample subset_n0 nulls and subset_n1 alternatives
  subset_idx = c(
    sample(idx_1, n1),
    sample(idx_0, n0)
  )
  subset_fold_change = data$fold_change[subset_idx,]
  subset_expression = data$expression[subset_idx,]

  subset_t_statistics = rowttests(
    subset_expression,
    colnames(subset_expression)
  )
  subset_true_Fdr = get_true_Fdr(
    subset_t_statistics$statistic,
    subset_fold_change$DE
  )

  subset = list(
    expression = subset_expression,
    statistics = subset_t_statistics$statistic,
    Fdr = subset_true_Fdr,
    fold_change = subset_fold_change
  )
  return(subset)
}

get_all_metrics <- function(
    method, estimated_fdr,
    estimated_pi0, data,
    cutoffs = c(0.2)
) {
  all_classification_metrics <- lapply(
    names(cutoffs),
    function(cutoff_name) {
      cutoff <- cutoffs[cutoff_name]
      metrics_cutoff <- classification_metrics(
        method,
        estimated_fdr,
        estimated_pi0,
        data$statistics,
        data$fold_change$DE,
        cutoff = cutoff
      )
      names(metrics_cutoff) <- paste0(
        names(metrics_cutoff), '_', cutoff_name
      )
      metrics_cutoff
    }
  ) %>%
    do.call(c, .)

  fdr_metrics <- method_metrics(
    method,
    estimated_fdr,
    data$statistics,
    data$fold_change$DE,
    data$Fdr
  ) %>% as.list()

  calib_metrics <- get_calib_metrics(
    estimated_fdr,
    data$fold_change$DE
  )

  fdr_metrics$pi0hat <- estimated_pi0

  c(all_classification_metrics, fdr_metrics, calib_metrics)
}


rowttests <- function(
    matrix,
    condition,
    var.equal = TRUE
) {
  condition.first = unique(condition)[1]
  data.frame(
    statistic = sapply(1:nrow(matrix), function(i) {
      row = matrix[i,]
      row1 = row[condition == condition.first]
      row2 = row[condition != condition.first]
      row_test = t.test(row1, row2, var.equal=var.equal)
      unname(row_test$statistic)
    })
  )
}

p_from_t <- function(test_statistics, df = NULL) {

  if (is.null(df)) {
    # assume standard normal
    one_sided <- unlist(lapply(test_statistics, function(z) {
      pnorm(-1*abs(z))
    }))
  } else {
    # assume t_df
    one_sided <-unlist(lapply(test_statistics, function(z) {
      pt(-1*abs(z), df = df)
    }))
  }

  return (2*one_sided)

}

get_Fdr_condition <- function(test_statistics, t, how = "symmetric") {
  if (!(how %in% c("quantile", "symmetric"))) {
    stop("how must be one of c(\"quantile\", \"symmetric\")")
  }

  if (how == "symmetric") {
    condition = abs(test_statistics) >= t
  } else if (how == "quantile") {
    if (t > median(test_statistics)) {
      percentile = mean(test_statistics > t)
      limits = quantile(test_statistics, c(percentile, 1-percentile))
    } else {
      percentile = mean(test_statistics < t)
      limits = quantile(test_statistics, c(percentile, 1-percentile))
    }
    condition = (
      test_statistics <= limits[1] |
        test_statistics >= limits[2]
    )
  }

  return(condition)
}

get_true_Fdr <- function(test_statistics, truth, how = "symmetric")  {
  out <- rep(NA, length(test_statistics))
  for (i in seq_len(length(test_statistics))) {
    t = test_statistics[i]
    condition = get_Fdr_condition(test_statistics, t, how = how)

    out[i] <- mean(1 - truth[condition])
  }
  return(out)
}

Fdr_from_fdr <- function(fdr, test_statistics, how = "symmetric") {
  Fdr = vector(length = length(test_statistics))
  for (i in seq_len(length(test_statistics))) {
    t = test_statistics[i]
    condition = get_Fdr_condition(test_statistics, t, how = how)

    Fdr[i] = mean(fdr[condition])
  }

  return(Fdr)
}

method_metrics <- function(
    method_name,
    estimated_fdr,
    test_statistics,
    hypothesis_labels,
    true_Fdr,
    estimated_Fdr = NULL,
    how_Fdr = "symmetric"
) {
  if (is.null(estimated_Fdr)) {
    estimated_Fdr = Fdr_from_fdr(estimated_fdr, test_statistics)
  }
  c(
    'method' = method_name,
    'roc' = get_roc(
      fdr = estimated_fdr,
      truth = hypothesis_labels
    ),
    'pr' = get_prauc(
      fdr = estimated_fdr,
      truth = hypothesis_labels
    ),
    'brier' = get_brier(
      fdr = estimated_fdr,
      truth = hypothesis_labels
    ),
    'Fdr MSE' = get_MSE(
      estimate = estimated_Fdr,
      true = true_Fdr
    )
  )
}

classification_metrics <- function(
    method, fdr, pi0, test_statistics, truth, cutoff = 0.2
) {
  pred = rep(1, length(test_statistics))
  pred[fdr > cutoff] = 0

  TP = sum(pred == 1 & truth == 1)
  FP = sum(pred == 1 & truth == 0)
  TN = sum(pred == 0 & truth == 0)
  FN = sum(pred == 0 & truth == 1)

  accuracy = (TP + TN)/(TP + TN + FP + FN)
  precision = TP/(TP+FP)
  sensitivity = TP/(TP+FN)
  specificity = TN/(TN + FP)
  FDR = FP/(FP + TP)
  f1 = (2 * precision * sensitivity) / (precision+sensitivity)
  prop_T = mean(pred == 1)

  return(list(
    'method' = method,
    'cutoff' = cutoff,
    'global_FDR' = FDR,
    'sensitivity' = sensitivity,
    'specificity' = specificity,
    'prop_pred_T' = prop_T,

    'TP' = TP,
    'FP' = FP,
    'TN' = TN,
    'FN' = FN,
    'accuracy' = accuracy,
    'precision' = precision,
    'f1' = f1
  ))
}

get_calib_metrics <- function(fdr, truth) {
    upper =  0.9999999
    lower = 0.0000001
    
    pr_1 = 1 - fdr # predicted probability that truth = 1
    
    pr_1[pr_1 > upper] = upper
    pr_1[pr_1 < lower] = lower
    log_pr_odds = log(pr_1/(1-pr_1)) # predicted log odds that truth = 1
    
    # "The log odds of predictions are used as the predictor of the 0/1 outcome"
    # in Ewout W. Steyerberg* and Yvonne Vergouwe 2014
    truth = as.numeric(truth)
    fit <- glm(
        truth~log_pr_odds,
        family = binomial(link='logit')
    )

    list(
        'calib_intercept' = unname(fit$coefficients['(Intercept)']),
        'calib_slope' = unname(fit$coefficients['log_pr_odds'])
    )
}

get_calib_dat <- function(
    fdrs, truth
) {
  calib_dat <- data.frame(
    fdrs
  )
  calib_dat$true = truth

  calib_dat <- calib_dat %>%
    pivot_longer(1:(ncol(calib_dat)-1), names_to = 'method', values_to = 'fdr') %>%
    mutate(
      fdr_group = cut(
        fdr,
        c(0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,1)
      )
    ) %>%
    group_by(method, fdr_group) %>%
    summarize(
      prop_null = 1 - mean(true),
      n = length(true),
      mean_fdr = mean(fdr),
      .groups = "keep"
    ) %>%
    mutate(
      fdr_group_number = case_when(
        fdr_group == "(0,0.01]" ~ 0.005,
        fdr_group == "(0.01,0.1]" ~ 0.05,
        fdr_group == "(0.1,0.2]" ~ 0.15,
        fdr_group == "(0.2,0.3]" ~ 0.25,
        fdr_group == "(0.3,0.4]" ~ 0.35,
        fdr_group == "(0.4,0.5]" ~ 0.45,
        fdr_group == "(0.5,0.6]" ~ 0.55,
        fdr_group == "(0.6,0.7]" ~ 0.65,
        fdr_group == "(0.7,0.8]" ~ 0.75,
        fdr_group == "(0.8,0.9]" ~ 0.85,
        fdr_group == "(0.9,0.99]" ~ 0.95,
        fdr_group == "(0.99,1]" ~ 0.995
      )
    )

  # add 95% CI based on binomial model
  z = qnorm(0.975)
  calib_dat <- calib_dat %>%
    mutate(
      upper = prop_null + z * sqrt(prop_null * (1 - prop_null) / n ),
      lower = prop_null - z * sqrt(prop_null * (1 - prop_null) / n )
    )
  return(calib_dat)
}

plot_calibration <- function(
    fdrs, truth, y_axis = "fdr_group_number"
) {
  calib_dat <- get_calib_dat(
    fdrs, truth
  )

  calib_dat <- calib_dat %>%
    mutate(
      fdr_group_number = ifelse(fdr_group == "(0.99,1]", 1, fdr_group_number),
      fdr_group_number = ifelse(fdr_group == "(0,0.01]", 0, fdr_group_number)
    ) %>%
    filter(!is.na(fdr_group))

  p <- calib_dat %>%
    ggplot(aes(x = !!rlang::sym(y_axis), y = prop_null, color = method)) +
    geom_segment(aes(x = 0, xend = 0.01, y = 0, yend = 0.01), lty = 2, col = "grey") +
    geom_segment(aes(x = 0.01, y = 0.01, xend = 0.99, yend = 0.99), lty = 2, col = 'grey') +
    geom_segment(aes(x = 0.99, xend = 1, y = 0.99, yend = 1), lty = 2, col = "grey") 
  
  if (y_axis == 'fdr_group_number') {
    p <- p +
      geom_point(size = 2, position = position_dodge(width = 0.05)) +
      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.05))
  } else {
    p <- p +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lower, ymax = upper))
  }
  p <- p +
    scale_color_manual(
      breaks = names(color_list),
      values = unlist(color_list)
    ) +
    scale_x_continuous(
      breaks = unique(calib_dat$fdr_group_number),
      labels = unique(calib_dat$fdr_group)
    ) +
    labs(
      x = "Binned Predicted Probability Null",
      y = "True Proportion Null",
      color = ""
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.12, 0.87),
      legend.title = element_text(size = 0.1),
      legend.box.background = element_rect(color= 'black'),
      axis.text.x = element_text(angle = 30,  hjust = 1)
    )

  return(p)
}

library(PRROC)
library(pROC)
library(dplyr)

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

get_MSE <- function(estimate, true) {
  mean((true - estimate)^2)
}

calibrate <- function(lfdr, true, plot = F) {
  yhat <- 1-lfdr
  calib_dat <- data.frame(
    true = as.numeric(as.logical(true)), prob = yhat, i = 1
  ) %>%
    arrange(-prob) %>%
    mutate(
      bin = as.numeric(case_when(
        prob < 0.1 ~ 1,
        prob < 0.2 ~ 2,
        prob < 0.3 ~ 3,
        prob < 0.4 ~ 4,
        prob < 0.5 ~ 5,
        prob < 0.6 ~ 6,
        prob < 0.7 ~ 7,
        prob < 0.8 ~ 8,
        prob < 0.9 ~ 9,
        TRUE ~ 10
      ))
    ) %>%
    group_by(bin) %>%
    summarize(
      proptrue = mean(true),
      avgprob = mean(prob),
      n = sum(i)
    )

  for (bin in 1:10) {
    if (!(bin %in% calib_dat$bin)) {
      calib_dat <- rbind(calib_dat, c(
        'bin' = bin, 'proptrue' = NA, 'avgprob' = NA, n = 0
      )) %>%
        arrange(bin)
    }
  }

  if (plot) {
    g <- calib_dat %>%
      ggplot(
        aes(x = avgprob, y = proptrue)
      ) +
      geom_point() +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, color = 'red') +
      ylim(0, 1) +
      scale_x_continuous(
        breaks = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95),
        labels = c('[0,0.1)', '[0.1,0.2)', '[0.2, 0.3)',
                   '[0.3, 0.4)', '[0.4, 0.5)', '[0.5, 0.6)',
                   '[0.6, 0.7)', '[0.7, 0.8)', '[0.8, 0.9)', '[0.9, 1)')
      ) +
      xlab('P(alternative) = 1-fdr') +
      ylab('Proportion alternative')
    print(g)
  }

  return(calib_dat)
}

calibrate_wide <- function(lfdr, true, plot = T) {
  calib_dat <- calibrate(lfdr, true, plot)
  calib_wide = c(calib_dat$proptrue, calib_dat$avgprob, calib_dat$n)
  names(calib_wide) = c(
    paste('proptrue', calib_dat$bin, sep = ''),
    paste('avgprob', calib_dat$bin, sep = ''),
    paste('n', calib_dat$bin, sep = '')
  )
  return(data.frame(t(calib_wide)))
}
