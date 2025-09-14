color_list = list(
  "fdrSAFE" = "#D55E00",
  "locfdr" = "#56B4E9",
  "fdrtool" = "#009E73",
  "qvalue" = "#0072B2",
  "fdrSAFE_selection-only" = "#E69F00",
  "fdrSAFE_ensemble-only" = "#CC79A7",
  "fdrSAFE_ensemble-all" = "#F0E442",
  "oracle" = "#000000",
  "oracle ensemble" = "#999999"
)

color_labels = list(
  "fdrSAFE" = "fdrSAFE",
  "locfdr" = "locfdr",
  "fdrtool" = "fdrtool",
  "qvalue" = "qvalue",
  "fdrSAFE_selection-only" = bquote(fdrSAFE[selection-only]),
  "fdrSAFE_ensemble-only" = bquote(fdrSAFE[ensemble-only]),
  "fdrSAFE_ensemble-all" = bquote(fdrSAFE[ensemble-all]),
  "oracle" = "oracle",
  "oracle ensemble" = "oracle ensemble"
)

load_res <- function(results_dir, study) {
    fdrSAFE_metrics <- read.csv(file.path(results_dir, "fdrSAFE_metrics.csv"))
    fdrSAFE_metrics$method = 'fdrSAFE'
    fdrSAFE_metrics$i = 1:nrow(fdrSAFE_metrics)

    grid_metrics <- read.csv(file.path(results_dir, "grid_metrics.csv"))
    grid_metrics$method = 'fdrSAFE_selection-only'
    grid_metrics$i = 1:nrow(grid_metrics)

    ensemble_metrics <- read.csv(file.path(results_dir, "ensemble_metrics.csv"))
    ensemble_metrics$method = 'fdrSAFE_ensemble-only'
    ensemble_metrics$i = 1:nrow(ensemble_metrics)

    ensemble_all_metrics <- read.csv(file.path(results_dir, "ensemble_all_metrics.csv"))
    ensemble_all_metrics$method = 'fdrSAFE_ensemble-all'
    ensemble_all_metrics$i = 1:nrow(ensemble_all_metrics)

    oracle_top_metrics <- read.csv(file.path(results_dir, "oracle_top_metrics.csv"))
    oracle_top_metrics$method = 'oracle'
    oracle_top_metrics$i = 1:nrow(oracle_top_metrics)
    oracle_top_metrics <- oracle_top_metrics[, colnames(fdrSAFE_metrics)]

    oracle_ensemble_metrics <- read.csv(file.path(results_dir, "oracle_ensemble_metrics.csv"))
    oracle_ensemble_metrics$method = 'oracle ensemble'
    oracle_ensemble_metrics$i = 1:nrow(oracle_ensemble_metrics)
    oracle_ensemble_metrics <- oracle_ensemble_metrics[, colnames(fdrSAFE_metrics)]

    qvalue_metrics <- read.csv(file.path(results_dir, "qvalue_metrics.csv"))
    qvalue_metrics$method = 'qvalue'
    qvalue_metrics$i = 1:nrow(qvalue_metrics)

    locfdr_metrics <- read.csv(file.path(results_dir, "locfdr_metrics.csv"))
    locfdr_metrics$method = 'locfdr'
    locfdr_metrics$i = 1:nrow(locfdr_metrics)

    fdrtool_metrics <- read.csv(file.path(results_dir, "fdrtool_metrics.csv"))
    fdrtool_metrics$method = 'fdrtool'
    fdrtool_metrics$i = 1:nrow(fdrtool_metrics)

    metrics = do.call(rbind, list(
        fdrSAFE_metrics, 
        grid_metrics,
        ensemble_metrics,
        ensemble_all_metrics,
        oracle_top_metrics,
        oracle_ensemble_metrics,
        qvalue_metrics,
        locfdr_metrics, 
        fdrtool_metrics
    ))
    metrics$study = study
    if (!('fdrerror' %in% colnames(metrics))) {
        metrics$fdrerror = metrics$Fdrerror
        metrics$fdrerror_topq = metrics$Fdrerror_topq
    }

    pi0 <- read.csv(file.path(results_dir, "pi0_estimates.csv"))
    pi0 <- pi0 %>%
        mutate(i = 1:nrow(pi0)) %>%
        pivot_longer(1:ncol(pi0), names_to = 'method', values_to = 'pi0') %>%
        mutate(method = case_when(
            method == 'grid' ~ 'fdrSAFE_selection-only',
            method == 'ensemble' ~ 'fdrSAFE_ensemble-only',
            method == 'ensemble_all' ~ 'fdrSAFE_ensemble-all',
            method == 'oracle_ensemble' ~ 'oracle ensemble',
            TRUE ~ method
        ))

    metrics <- pi0 %>%
        merge(metrics, by = c('i', 'method'))

    return(metrics)
}