library(kcde)
library(pdtmvn)
library(plyr)
library(dplyr)
library(lubridate)

options(error = recover)

source("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/code/sim-densities-sim-study-discretized-Duong-Hazelton.R")


all_data_sets <- c("sim")
all_bw_parameterizations <- c("diagonal", "full")
all_sim_n <- c("100", "1000")
all_sim_families <- "multivariate-2d-discretized"
all_sim_run_inds <- seq(from = 1, to = 500)



## Test values for debugging
#all_data_sets <- c("ili_national")
#all_prediction_horizons <- 1L
#all_max_lags <- 1L
#all_filtering_values <- c(TRUE)
#all_seasonality_values <- c(TRUE)
#all_bw_parameterizations <- c("full")

## Allocate data frame to store results
num_rows <- length(all_sim_n) *
    length(all_sim_families) *
    length(all_sim_run_inds) *
    length(all_bw_parameterizations)

results <- data.frame(
    bw_parameterization = rep(NA_character_, num_rows),
    sim_n = rep(NA_character_, num_rows),
    sim_family = rep(NA_character_, num_rows),
    sim_run_inds = rep(NA_character_, num_rows),
    model = "kcde",
    KL_div = rep(NA_real_, num_rows),
    ISE = rep(NA_real_, num_rows),
    stringsAsFactors = FALSE
)

results_row_ind <- 1L

n_eval <- 10^5
for(sim_family in all_sim_families) {
    ## Set path where fit object is stored
    results_path <- file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/sim-study",
        sim_family,
        "estimation-results")
    
    ## Draw a sample from true density that was being estimated
    sample_true_density <- sim_from_pdtmvn_mixt(n = n_eval, sim_family = sim_family)
    true_log_density_sample_true_density <- d_pdtmvn_mixt_conditional(X = sample_true_density,
        sim_family = sim_family,
        conditional = TRUE,
        log = TRUE)
    
    for(sim_run_ind in all_sim_run_inds) {
        for(sim_n in all_sim_n) {
            for(bw_parameterization in all_bw_parameterizations) {
                ## Set values describing case in results
                results$sim_family[results_row_ind] <- sim_family
                results$sim_run_ind[results_row_ind] <- sim_run_ind
                results$sim_n[results_row_ind] <- sim_n
                results$bw_parameterization[results_row_ind] <-
                    bw_parameterization
                
                ## Load kcde_fit object.  Estimation was performed previously.
                case_descriptor <- paste0(
                    "sim",
                    "_", sim_run_ind,
                    "-prediction_horizon_", 0L,
                    "-max_lag_", 0L,
                    "-max_seasonal_lag_", 0L,
                    "-filtering_", "FALSE",
                    "-differencing_", "FALSE",
                    "-seasonality_", "FALSE",
                    "-bw_parameterization_", bw_parameterization,
                    "-sim_n_", sim_n,
                    "-sim_ind_", sim_run_ind)
                
                
                kcde_fit <- readRDS(file.path(results_path,
                        paste0("kcde_fit-", case_descriptor, ".rds")))
                
                ## Estimated density evaluated at the sample drawn from the true density above
                unique_sample_true_density <- sample_true_density %>%
                    as.data.frame %>%
                    distinct
                
                unique_sample_true_density$est_log_density_unique_sample_true_density <- apply(unique_sample_true_density, 1, function(data_row) {
                        kcde_predict(
                            kcde_fit = kcde_fit,
                            prediction_data =
                                matrix(data_row[seq_len(length(data_row) - 1)],
                                    nrow = 1,
                                    dimnames = list(NULL, names(data_row[seq_len(length(data_row) - 1)]))),
                            leading_rows_to_drop = 0L,
                            trailing_rows_to_drop = 0L,
                            additional_training_rows_to_drop = NULL,
                            prediction_type = "distribution",
                            prediction_test_lead_obs = matrix(data_row[length(data_row)],
                                nrow = 1,
                                dimnames = list(NULL, paste0(names(data_row[length(data_row)]), "_horizon0"))),
                            log = TRUE
                        )
                    })
                
                est_log_density_sample_true_density <- sample_true_density %>%
                    as.data.frame %>%
                    left_join(unique_sample_true_density, by = colnames(sample_true_density))
                
                est_log_density_sample_true_density <- est_log_density_sample_true_density$est_log_density_unique_sample_true_density
                
                
                ## Estimate K-L Divergence
                results$KL_div[results_row_ind] <- mean(true_log_density_sample_true_density - est_log_density_sample_true_density)
                
                ## Estimate Hellinger Distance
                results$Hellinger_dist[results_row_ind] <-
                    sqrt(1 - mean(exp(0.5 * (est_log_density_sample_true_density - true_log_density_sample_true_density))))
                
                results_row_ind <- results_row_ind + 1
            } # bw_parameterization
        } # sim_n
    } # sim_run_ind
    
    saveRDS(sample_true_density, file = paste0(
            "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/sim-study/sample_true_density_for_kcde-predictions_",
            sim_family,
            ".rds"))
} # sim_family

## Save results
saveRDS(results, file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/sim-study/kcde-predictions.rds"))
