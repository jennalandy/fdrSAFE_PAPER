library(fdrSAFE)
library(fdrtool)
library(locfdr)
library(qvalue)
library(stringr)

source("test_eval.R")

# print(fdrSAFE)
run_test <- function(
    data_generating_fn, 
    N, 
    path,
    ensemble_size = 10, 
    n_synthetic = 10,
    overwrite = FALSE, 
    sim_n = NULL,
    parallel = FALSE,
    locfdr_grid = NULL,
    fdrtool_grid = NULL,
    qvalue_grid = NULL,
    df = NULL,
    to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)},
    type = 'symmetric',
    focus_metric = 'fdrerror'
) {
  redo_pval_function = FALSE
  how = ""
  if (class(to_pval_function) == 'character') {
    how = to_pval_function
    redo_pval_function = TRUE
  }

  given_locfdr_grid = locfdr_grid
  given_fdrtool_grid = fdrtool_grid
  given_qvalue_grid = qvalue_grid
  
  if (!endsWith(path, '/')) {
    path = paste(path, '/', sep = '')
  }
  
  if (!dir.exists(path)) {
    dir.create(path)
  } else if (overwrite) {
    cat('\nDirectory', path, 'already exists')
    cat('\nOverwritting because given overwrite = TRUE')
  } else {
    cat('\nDirectory', path, 'already exists')
    cat('\nEnding tests because given overwrite = FALSE')
    cat('\nRerun with overwrite = TRUE or choose a different directory\n\n')
    return()
  }
  
  sink(file = paste(path, 'tests.log', sep = ''))
  
  start = Sys.time()
  cat(paste("Starting Test at", start))
  rep = 1
  while (rep <= N) {
    cat(paste("\n\nRep", rep, '/', N))
    tryCatch({
      dat <- data_generating_fn()
      if (redo_pval_function) {
        if (how == "empirical_null_var") {
          empirical_null_var = var(dat$zz[dat$truth == 0])
          print(empirical_null_var)
          to_pval_function = function(test_statistics) {
            sapply(test_statistics, function(t) {
              2*(1-pnorm(abs(t), 0, sqrt(empirical_null_var)))
            })
          }
        } else if (how == "empirical_var") {
          empirical_var = var(dat$zz)
          print(empirical_var)
          to_pval_function = function(test_statistics) {
            sapply(test_statistics, function(t) {
              2*(1-pnorm(abs(t), 0, sqrt(empirical_var)))
            })
          }
        } else if (how == 'perm_var') {
          to_pval_function = function(test_statistics) {
            sapply(test_statistics, function(t) {
              2*(1-pnorm(abs(t), 0, sqrt(dat$perm_var)))
            })
          }
        }
      }
      
      if (is.null(sim_n)) {
        sim_n = length(dat$zz)
      }
      
      if (is.null(given_locfdr_grid) & 
          is.null(given_fdrtool_grid) & 
          is.null(given_qvalue_grid)) {
        locfdr_grid <- build_locfdr_grid(dat$zz)
        fdrtool_grid <- build_fdrtool_grid(dat$zz)
        qvalue_grid <- build_qvalue_grid(dat$zz, to_pval_function = to_pval_function)
      } else {
        locfdr_grid <- given_locfdr_grid
        fdrtool_grid <- given_fdrtool_grid
        qvalue_grid <- given_qvalue_grid
      }
      
      # make sure defaults are included
      qvalue_grid <- rbind(
        qvalue_grid,
        c('probit', 1.5, 'smoother', FALSE)
      ) %>%
        unique()
      qvalue_grid <- data.frame(qvalue_grid)
      qvalue_grid$adj <- as.numeric(qvalue_grid$adj)
      
      fdrtool_grid <- rbind(
        fdrtool_grid,
        c("fndr", 0.75)
      ) %>%
        unique()
      
      locfdr_grid <- rbind(
        locfdr_grid,
        c(0, 1/4, 1, 0)
      ) %>%
        unique()
      
      
      locfdr_grid <- reduce_locfdr_grid(dat$zz, locfdr_grid)
      fdrtool_grid <- reduce_fdrtool_grid(dat$zz, fdrtool_grid)
      qvalue_grid <- reduce_qvalue_grid(dat$zz, qvalue_grid = qvalue_grid, to_pval_function = to_pval_function)
      
      locfdr_default_row <- locfdr_grid %>%
        mutate(row = 1:nrow(locfdr_grid)) %>%
        filter(
          pct==0,
          pct0==1/4,
          nulltype == 1,
          type == 0
        ) %>%
        pull(row)
      fdrtool_default_row <- fdrtool_grid %>%
        mutate(row = 1:nrow(fdrtool_grid)) %>%
        filter(
          cutoff.method == 'fndr',
          pct0 == 0.75
        ) %>%
        pull(row)
      qvalue_default_row <- qvalue_grid %>%
        mutate(row = 1:nrow(qvalue_grid)) %>%
        filter(
          transf == 'probit',
          adj == 1.5,
          pi0.method == 'smoother',
          as.logical(smooth.log.pi0) == FALSE
        ) %>%
        pull(row)
      
      grid_size = nrow_null0(locfdr_grid) +
        nrow_null0(fdrtool_grid) +
        nrow_null0(qvalue_grid)
      cat(paste('\n\tgrid grid_size =', grid_size))
      
      if (ensemble_size > grid_size) {
        cat('\n\tensemble size > a grid size, skipping')
        stop()
      }
      
      topq <- abs(dat$zz) > quantile(abs(dat$zz), 0.75)
      
      cat("\n\tsetup done")
      
      # ------- ORACLE ------- 
      oracle_grid <- get_oracle_grid(
        dat = dat, 
        to_pval_function = to_pval_function, 
        grids = list(
          'locfdr' = locfdr_grid,
          'fdrtool' = fdrtool_grid,
          'qvalue' = qvalue_grid
        ),
        topq = topq
      )
            
      # ------- RUN ALL METHODS ------- 
      
      # default runs of other methods
      locfdr_run <- locfdr::locfdr(dat$zz, plot = F)
      locfdr_run$Fdr <- Fdr_from_fdr(locfdr_run$fdr, dat$zz)
      
      fdrtool_run <- fdrtool::fdrtool(dat$zz, plot = F, verbose = F)
      fdrtool_run$Fdr <- Fdr_from_fdr(fdrtool_run$lfdr, dat$zz)
      
      qvalue_run <- NULL
      tryCatch({
        qvalue_run <- qvalue::qvalue(
          to_pval_function(dat$zz), plot = 0, verbose = F
        )
      }, error = function(e) {})
      if (is.null(qvalue_run)) {
        qvalue_run <- qvalue::qvalue(
          to_pval_function(dat$zz), plot = 0, lambda = 0, verbose = F
        )
      }
      qvalue_run$Fdr <- Fdr_from_fdr(qvalue_run$lfdr, dat$zz)
      
      cat("\n\t\tdefaults run")
      oracle_focus_metric = focus_metric
      if (!(oracle_focus_metric %in% colnames(oracle_grid))) {
        oracle_focus_metric = str_replace(oracle_focus_metric, 'fdrerror','Fdrerror')
      }
      
      oracle_top_ensemble_idx <- head(oracle_grid[order(oracle_grid[,oracle_focus_metric]),], ensemble_size)
      oracle_locfdr_grid <- locfdr_grid[
        oracle_top_ensemble_idx[oracle_top_ensemble_idx$method == 'locfdr', 'row'],
      ]
      oracle_fdrtool_grid <- fdrtool_grid[
        oracle_top_ensemble_idx[oracle_top_ensemble_idx$method == 'fdrtool', 'row'],
      ]
      oracle_qvalue_grid <- qvalue_grid[
        oracle_top_ensemble_idx[oracle_top_ensemble_idx$method == 'qvalue', 'row'],
      ]
    
      oracle_ensemble_run <- fdrSAFE(
        dat$zz, n_synthetic = 0, n_workers = 5,
        ensemble_size = ensemble_size, verbose = F,
        locfdr_grid = oracle_locfdr_grid,
        fdrtool_grid = oracle_fdrtool_grid,
        qvalue_grid = oracle_qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df,
        to_pval_function = to_pval_function,
        type = type
      )
      
      cat("\n\t\toracle ensemble run")
      
      ensemble_run <- fdrSAFE(
        dat$zz, n_synthetic = 0, n_workers = 5,
        ensemble_size = ensemble_size, verbose = F,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df,
        to_pval_function = to_pval_function,
        type = type
      )
      cat("\n\t\tensemble run")
      
      ensemble_all_run <- fdrSAFE(
        dat$zz, n_synthetic = 0, n_workers = 5,
        ensemble_size = grid_size, verbose = F,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df,
        to_pval_function = to_pval_function,
        type = type
      )
      
      cat("\n\t\tensemble all run")
      
      grid_run <- fdrSAFE(
        dat$zz, n_synthetic = n_synthetic, 
        ensemble_size = 1,
        verbose = F, n_workers = 5,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df,
        to_pval_function = to_pval_function,
        type = type
      )
      
      cat("\n\t\tgrid run")
      
      fdrSAFE_run <- fdrSAFE(
        dat$zz, n_synthetic = n_synthetic,
        ensemble_size = ensemble_size,
        verbose = F, 
        n_workers = 5,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df,
        to_pval_function = to_pval_function,
        type = type
      )
      
      cat("\n\t\tfdrSAFE run")
      
      cat("\n\tall methods run")
      
      # ------- COMPUTE METRICS ------- 
      
      oracle_top_metrics_i <- head(oracle_grid[order(oracle_grid[,oracle_focus_metric]),], 1)
      
      oracle_med_metrics_i <- oracle_grid[order(oracle_grid[,oracle_focus_metric]),][round(nrow(oracle_grid)/2),]
      
      oracle_worst_metrics_i <- head(oracle_grid[order(oracle_grid[,oracle_focus_metric], decreasing = TRUE),], 1)
      
      
      generate_params_i <- fdrSAFE_run$working_model$parameters

      pi0_estimates_i <- list(
        "locfdr" = locfdr_run$fp0["mlest", "p0"],
        "fdrtool" = fdrtool_run$param[1, 'eta0'],
        "qvalue" = qvalue_run$pi0,
        "oracle_top" = oracle_top_metrics_i$pi0,
        "oracle_med" = oracle_med_metrics_i$pi0,
        "oracle_ensemble" = oracle_ensemble_run$pi0,
        "ensemble" = ensemble_run$pi0,
        "ensemble_all" = ensemble_all_run$pi0,
        "grid" = grid_run$pi0,
        "fdrSAFE" = fdrSAFE_run$pi0
      )
      
      pi0_var_i <- list(
        "oracle_ensemble" = oracle_ensemble_run$pi0_var,
        "ensemble" = ensemble_run$pi0_var,
        "ensemble_all" = ensemble_all_run$pi0_var,
        "fdrSAFE" = fdrSAFE_run$pi0_var
      )
      
      fdr_var_i <- list(
        "oracle_ensemble" = median(oracle_ensemble_run$fdr_var),
        "ensemble" = median(ensemble_run$fdr_var),
        "ensemble_all" = median(ensemble_all_run$fdr_var),
        "fdrSAFE" = median(fdrSAFE_run$fdr_var)
      )
      
      locfdr_cutoff <- quantile(locfdr_run$fdr, 1 - min(pi0_estimates_i$locfdr, 1))
      locfdr_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = locfdr_run$fdr, 
        Fdr = locfdr_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = locfdr_cutoff
      )
      
      fdrtool_cutoff <- quantile(fdrtool_run$lfdr, 1 - min(pi0_estimates_i$fdrtool, 1))
      fdrtool_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = fdrtool_run$lfdr, 
        Fdr = fdrtool_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = fdrtool_cutoff
      )
      
      qvalue_cutoff <- quantile(qvalue_run$lfdr, 1 - min(pi0_estimates_i$qvalue, 1))
      qvalue_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = qvalue_run$lfdr, 
        Fdr = qvalue_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = qvalue_cutoff
      )
      
      oracle_ensemble_cutoff <- quantile(oracle_ensemble_run$fdr, 1 - min(pi0_estimates_i$oracle_ensemble))
      oracle_ensemble_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = oracle_ensemble_run$fdr, 
        Fdr = oracle_ensemble_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = oracle_ensemble_cutoff
      )
      
      ensemble_cutoff <- quantile(ensemble_run$fdr, 1 - min(pi0_estimates_i$ensemble))
      ensemble_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = ensemble_run$fdr, 
        Fdr = ensemble_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = ensemble_cutoff
      )
      
      ensemble_all_cutoff <- quantile(ensemble_all_run$fdr, 1 - min(pi0_estimates_i$ensemble_all))
      ensemble_all_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = ensemble_all_run$fdr, 
        Fdr = ensemble_all_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = ensemble_all_cutoff
      )
      
      grid_cutoff <- quantile(grid_run$fdr, 1 - min(pi0_estimates_i$grid))
      grid_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = grid_run$fdr, 
        Fdr = grid_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = grid_cutoff
      )
      
      fdrSAFE_cutoff <- quantile(fdrSAFE_run$fdr, 1 - min(pi0_estimates_i$fdrSAFE))
      fdrSAFE_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = fdrSAFE_run$fdr, 
        Fdr = fdrSAFE_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = fdrSAFE_cutoff
      )
      
      fdrSAFE_run$top_grid$method = factor(fdrSAFE_run$top_grid$method, levels = c('qvalue','locfdr','fdrtool'))
      fdrSAFE_package_counts_i = table(fdrSAFE_run$top_grid$method)
      fdrSAFE_package_counts_i$i = rep
    
      final_grid <- fdrSAFE_run$all_grids %>%
        # select(-sim) %>%
        group_by(method, row) %>%
        summarize_all(mean)
      
      names(final_grid)[names(final_grid) == oracle_focus_metric] = paste0(oracle_focus_metric, "_est")
      oracle_final_grid <- merge(
        final_grid[,c('method','row', paste0(oracle_focus_metric, "_est"))],
        oracle_grid[,c('method','row',oracle_focus_metric)]
      )
      
      objective_bias = mean((oracle_final_grid[,paste0(oracle_focus_metric, "_est")] - oracle_final_grid[,oracle_focus_metric]), na.rm = TRUE)
      objective_MSE = mean((oracle_final_grid[,paste0(oracle_focus_metric, "_est")] - oracle_final_grid[,oracle_focus_metric])**2, na.rm = TRUE)
      objective_var = var((oracle_final_grid[,paste0(oracle_focus_metric, "_est")] - oracle_final_grid[,oracle_focus_metric]), na.rm = TRUE)
      objective_cor = cor(oracle_final_grid[,paste0(oracle_focus_metric, "_est")][!is.na(oracle_final_grid[,paste0(oracle_focus_metric, "_est")])], 
                          oracle_final_grid[,oracle_focus_metric][!is.na(oracle_final_grid[,paste0(oracle_focus_metric, "_est")])])
      
      oracle_final_grid <- oracle_final_grid[order(oracle_final_grid[,oracle_focus_metric]),]
      oracle_final_grid$oracle_rank <- 1:nrow(oracle_final_grid)
      oracle_final_grid <- oracle_final_grid[order(oracle_final_grid[,paste0(oracle_focus_metric, "_est")]),]
      chosen_oracle_rank_summary <- summary(head(oracle_final_grid, ensemble_size)$oracle_rank)
      
      names(chosen_oracle_rank_summary) = paste0('rank_', names(chosen_oracle_rank_summary))
      
      oracle_metrics_i <- c(list(
        'objective_bias' = objective_bias,
        'objective_MSE' = objective_MSE,
        'objective_var' = objective_var,
        'objective_cor' = objective_cor
      ), chosen_oracle_rank_summary)
      
      # STORE METRICS
      if (rep == 1) {
        generate_params <- data.frame(generate_params_i)
        pi0_estimates <- data.frame(pi0_estimates_i)
        pi0_var <- data.frame(pi0_var_i)
        fdr_var <- data.frame(fdr_var_i)
        
        locfdr_metrics <- data.frame(locfdr_metrics_i)
        fdrtool_metrics <- data.frame(fdrtool_metrics_i)
        qvalue_metrics <- data.frame(qvalue_metrics_i)
        
        ensemble_metrics <- data.frame(ensemble_metrics_i)
        ensemble_all_metrics <- data.frame(ensemble_all_metrics_i)
        grid_metrics <- data.frame(grid_metrics_i)
        fdrSAFE_metrics <- data.frame(fdrSAFE_metrics_i)
        
        fdrSAFE_package_counts <- data.frame(fdrSAFE_package_counts_i)
        
        oracle_worst_metrics = data.frame(oracle_worst_metrics_i)
        oracle_top_metrics = data.frame(oracle_top_metrics_i)
        oracle_med_metrics = data.frame(oracle_med_metrics_i)
        oracle_ensemble_metrics = data.frame(oracle_ensemble_metrics_i)
        oracle_metrics = data.frame(oracle_metrics_i)
      } else {
        generate_params <- rbind(
          generate_params,
          generate_params_i
        )
        pi0_estimates <- rbind(
          pi0_estimates, 
          pi0_estimates_i
        )
        pi0_var <- rbind(
          pi0_var,
          pi0_var_i
        )
        fdr_var <- rbind(
          fdr_var,
          fdr_var_i
        )
        
        locfdr_metrics <- rbind(
          locfdr_metrics, 
          locfdr_metrics_i
        )
        fdrtool_metrics <- rbind(
          fdrtool_metrics,
          fdrtool_metrics_i
        )
        qvalue_metrics <- rbind(
          qvalue_metrics,
          qvalue_metrics_i
        )
        
  
        ensemble_metrics <- rbind(
          ensemble_metrics,
          ensemble_metrics_i
        )
        ensemble_all_metrics <- rbind(
          ensemble_all_metrics,
          ensemble_all_metrics_i
        )
        grid_metrics <- rbind(
          grid_metrics,
          grid_metrics_i
        )
        fdrSAFE_metrics <- rbind(
          fdrSAFE_metrics,
          fdrSAFE_metrics_i
        )
        
        fdrSAFE_package_counts <- rbind(
          fdrSAFE_package_counts,
          fdrSAFE_package_counts_i
        )
        
        oracle_top_metrics <- rbind(
          oracle_top_metrics,
          data.frame(oracle_top_metrics_i)
        )
        oracle_med_metrics <- rbind(
          oracle_med_metrics, 
          data.frame(oracle_med_metrics_i)
        )
        oracle_worst_metrics <- rbind(
          oracle_worst_metrics,
          data.frame(oracle_worst_metrics_i)
        )
        oracle_ensemble_metrics <- rbind(
          oracle_ensemble_metrics,
          data.frame(oracle_ensemble_metrics_i)
        ) 
        oracle_metrics <- rbind(
          oracle_metrics,
          data.frame(oracle_metrics_i)
        )
      }
      
      cat("\n\tMETRICS computed")
      
      # ------- SAVE ALL DATAFRAMES ON EACH ITER (JUST IN CASE) ------- 
      write.csv(generate_params, paste0(path, 'generate_params', '.csv'), row.names = FALSE)
      write.csv(pi0_estimates, paste0(path, 'pi0_estimates', '.csv'), row.names = FALSE)
      write.csv(pi0_var, paste0(path, 'pi0_var', '.csv'), row.names = FALSE)
      write.csv(fdr_var, paste0(path, 'fdr_var', '.csv'), row.names = FALSE)
      
      
      write.csv(locfdr_metrics, paste0(path, 'locfdr_metrics', '.csv'), row.names = FALSE)
      write.csv(fdrtool_metrics, paste0(path, 'fdrtool_metrics', '.csv'), row.names = FALSE)
      write.csv(qvalue_metrics, paste0(path, 'qvalue_metrics', '.csv'), row.names = FALSE)

      
      write.csv(ensemble_metrics, paste0(path, 'ensemble_metrics', '.csv'), row.names = FALSE)
      write.csv(ensemble_all_metrics, paste0(path, 'ensemble_all_metrics', '.csv'), row.names = FALSE)
      write.csv(grid_metrics, paste0(path, 'grid_metrics', '.csv'), row.names = FALSE)
      write.csv(fdrSAFE_metrics, paste0(path, 'fdrSAFE_metrics', '.csv'), row.names = FALSE)
      
      write.csv(fdrSAFE_package_counts, paste0(path, 'fdrSAFE_package_counts', '.csv'), row.names = FALSE)
      
      write.csv(oracle_top_metrics, paste0(path, "oracle_top_metrics.csv"), row.names = FALSE)
      write.csv(oracle_med_metrics, paste0(path, "oracle_med_metrics.csv"), row.names = FALSE)
      write.csv(oracle_worst_metrics, paste0(path, "oracle_worst_metrics.csv"), row.names = FALSE)
      write.csv(oracle_ensemble_metrics, paste0(path, "oracle_ensemble_metrics.csv"), row.names = FALSE)
      write.csv(oracle_metrics, paste0(path, "oracle_metrics.csv"), row.names = FALSE)
      
      cat("\n\tdataframes saved")
      cat(paste("\nFinished", rep, '/', N, ': ', Sys.time()))
      
      rep = rep + 1
    },
    error = function(e) {
      cat(paste('\nError on', rep,'\n'))
      cat(paste("\n", e))
    })
  }
  
  now <- Sys.time()
  cat(paste('\n\nDone in', round(difftime(now, start, units = 'hours'), 2), 'hours'))
}



run_inc_n_test <- function(
    data_generating_fn, 
    N, 
    path,
    overwrite = FALSE,
    df = NULL,
    to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)},
    type = 'symmetric'
) {
  if (!endsWith(path, '/')) {
    path = paste(path, '/', sep = '')
  }
  
  if (!dir.exists(path)) {
    dir.create(path)
  } else if (overwrite) {
    cat('\nDirectory', path, 'already exists')
    cat('\nOverwritting because given overwrite = TRUE')
  } else {
    cat('\nDirectory', path, 'already exists')
    cat('\nEnding tests because given overwrite = FALSE')
    cat('\nRerun with overwrite = TRUE or choose a different directory\n\n')
    return()
  }
  
  sink(file = paste(path, 'tests.log', sep = ''))
  
  n_synthetic_values = c(0, 1, 5, 10, 20)
  ensemble_size_values = c(1, 5, 10, 20)
  
  start <- Sys.time()
  cat(paste("Starting test at", start))
  rep = 1
  while (rep <= N) {
    cat(paste('\n\nRep', rep, '/', N))
    
    dat <- data_generating_fn()
    
    locfdr_grid <- build_locfdr_grid(dat$zz)#, method = 'grid')#, lower_pi0 = lower_pi0)
    fdrtool_grid <- build_fdrtool_grid(dat$zz)#, method = 'grid')#, lower_pi0 = lower_pi0)
    qvalue_grid <- build_qvalue_grid(dat$zz)#, method = 'grid')#, lower_pi0 = lower_pi0)
    grid_size = nrow_null0(locfdr_grid) +
      nrow_null0(fdrtool_grid) +
      nrow_null0(qvalue_grid)
    cat(paste('\n\tGrid_size =', grid_size))
    
    if (grid_size < 20) {
      cat('\ngridsize too small (<20), skipping dataset')
      next
    }
    
    true_Fdr = dat$Fdr
    true_fdr = dat$fdr
    topq = abs(dat$zz) > quantile(abs(dat$zz), 0.75)
    
    for (n_synthetic in n_synthetic_values) {
      cat(paste('\n\tn_synthetic =', n_synthetic))
      
      for (ensemble_size in ensemble_size_values) {
        cat(paste('\n\t\tensemble_size =', ensemble_size))
        
        fdrSAFE_run <- fdrSAFE(
          test_statistics = dat$zz, 
          n_synthetic = n_synthetic,
          ensemble_size = ensemble_size,
          verbose = FALSE, 
          parallel = FALSE,
          locfdr_grid = locfdr_grid,
          fdrtool_grid = fdrtool_grid,
          qvalue_grid = qvalue_grid, 
          df = df,
          type = type,
          to_pval_function = to_pval_function
        )
        cat('\n\t\t\tfdrSAFE done')
        
        ensemble_run <-  fdrSAFE(
          dat$zz, 
          n_synthetic = 0,
          ensemble_size = ensemble_size, 
          verbose = FALSE,
          parallel = FALSE,
          locfdr_grid = locfdr_grid,
          fdrtool_grid = fdrtool_grid,
          qvalue_grid = qvalue_grid,
          df = df,
          type = type,
          to_pval_function = to_pval_function
        )
        cat('\n\t\t\tensemble done')
        
        fdrSAFE_cutoff <- quantile(fdrSAFE_run$fdr, 1 - min(fdrSAFE_run$pi0))
        fdrSAFE_metrics_i <- calc_metrics(
          test_statistics = dat$zz,
          fdr = fdrSAFE_run$fdr, 
          Fdr = fdrSAFE_run$Fdr, 
          truth = dat$truth, 
          true_Fdr = dat$true_Fdr,
          true_fdr = dat$true_fdr,
          topq = topq,
          cutoff = fdrSAFE_cutoff
        )
        fdrSAFE_metrics_i$n_synthetic = n_synthetic
        fdrSAFE_metrics_i$ensemble_size = ensemble_size
        
        ensemble_cutoff <- quantile(ensemble_run$fdr, 1 - min(ensemble_run$pi0))
        ensemble_metrics_i <- calc_metrics(
          test_statistics = dat$zz,
          fdr = ensemble_run$fdr, 
          Fdr = ensemble_run$Fdr, 
          truth = dat$truth, 
          true_Fdr = dat$true_Fdr,
          true_fdr = dat$true_fdr,
          topq = topq,
          cutoff = ensemble_cutoff
        )
        ensemble_metrics_i$n_synthetic = n_synthetic
        ensemble_metrics_i$ensemble_size = ensemble_size
        
        if (
          rep == 1 &
          n_synthetic == n_synthetic_values[1] & 
          ensemble_size == ensemble_size_values[1]
        ) {
          fdrSAFE_metrics <- data.frame(fdrSAFE_metrics_i)
          ensemble_metrics <- data.frame(ensemble_metrics_i)
        } else {
          fdrSAFE_metrics <- rbind(
            fdrSAFE_metrics,
            fdrSAFE_metrics_i
          )
          ensemble_metrics <- rbind(
            ensemble_metrics, 
            ensemble_metrics_i
          )
        }
        
        write.csv(
          fdrSAFE_metrics, 
          paste0(path, "fdrSAFE_metrics", ".csv")
        )
        write.csv(
          ensemble_metrics, 
          paste0(path, "ensemble_metrics", ".csv")
        )
        cat("\n\t\t\tdataframes saved")
      }
    }
    cat(paste("\n\t\tFinished", rep, '/', N, ': ', Sys.time()))
    rep = rep + 1
  }
  now <- Sys.time()
  cat(paste('\n\nDone in', round(difftime(now, start, units = 'hours'), 2), 'hours'))
  sink()
}