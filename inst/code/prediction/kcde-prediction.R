library(kcde)
library(pdtmvn)
library(plyr)
library(dplyr)
library(lubridate)

options(error = recover)

all_data_sets <- c("ili_national")
all_prediction_horizons <- seq_len(52)
all_max_lags <- 1L
all_filtering_values <- c(FALSE, TRUE)
all_seasonality_values <- c(FALSE, TRUE)
all_bw_parameterizations <- c("diagonal", "full")

## Test values for debugging
#all_data_sets <- c("ili_national")
#all_prediction_horizons <- 1L
#all_max_lags <- 1L
#all_filtering_values <- c(TRUE)
#all_seasonality_values <- c(TRUE)
#all_bw_parameterizations <- c("full")


for(data_set in all_data_sets) {
    ## Set path where fit object is stored
    results_path <- file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/results",
        data_set,
        "estimation-results")
    
    
    if(identical(data_set, "ili_national")) {
        ## Load data for nationally reported influenza like illness
        library(cdcfluview)
        
        usflu <- get_flu_data("national", "ilinet", years=1997:2014)
        data <- transmute(usflu,
            region.type = REGION.TYPE,
            region = REGION,
            year = YEAR,
            week = WEEK,
            weighted_ili = as.numeric(X..WEIGHTED.ILI))
        
        ## Subset data to do prediction using only data up through 2014
        data <- data[data$year <= 2014, , drop = FALSE]
        
        ## Row indices in data corresponding to times at which we want to make a prediction
        prediction_time_inds <- which(data$year %in% 2011:2014)
        
        ## Add time column.  This is used for calculating times to drop in cross-validation
        data$time <- ymd(paste(data$year, "01", "01", sep = "-"))
        week(data$time) <- data$week
        
        ## Add time_index column.  This is used for calculating the periodic kernel.
        ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
        ## The origin is arbitrary.
        data$time_index <- as.integer(data$time - ymd(paste("1970", "01", "01", sep = "-")))
        
        ## Set prediction target var
        prediction_target_var <- "weighted_ili"
    }
    
    
    ## Allocate data frame to store results for this data set
    num_rows <- length(all_prediction_horizons) *
        length(all_max_lags) *
        length(all_filtering_values) *
        length(all_seasonality_values) *
        length(all_bw_parameterizations) *
        length(prediction_time_inds)
    
    data_set_results <- data.frame(data_set = data_set,
        prediction_horizon = rep(NA_integer_, num_rows),
        max_lag = rep(NA_integer_, num_rows),
        filtering = rep(NA, num_rows),
        seasonality = rep(NA, num_rows),
        bw_parameterization = rep(NA_character_, num_rows),
        model = "kcde",
        prediction_time = rep(NA, num_rows),
        log_score = rep(NA_real_, num_rows),
        pt_pred = rep(NA_real_, num_rows),
        AE = rep(NA_real_, num_rows),
        interval_pred_lb_95 = rep(NA_real_, num_rows),
        interval_pred_ub_95 = rep(NA_real_, num_rows),
        interval_pred_lb_50 = rep(NA_real_, num_rows),
        interval_pred_ub_50 = rep(NA_real_, num_rows),
        stringsAsFactors = FALSE
    )
    class(data_set_results$prediction_time) <- class(data$time)
    
    results_row_ind <- 1L
    for(prediction_horizon in all_prediction_horizons) {
        for(max_lag in all_max_lags) {
            for(filtering in all_filtering_values) {
                for(seasonality in all_seasonality_values) {
                    for(bw_parameterization in all_bw_parameterizations) {
                        for(prediction_time_ind in prediction_time_inds) {
                            ## Set values describing case in data_set_results
                            data_set_results$prediction_horizon[results_row_ind] <-
                                prediction_horizon
                            data_set_results$max_lag[results_row_ind] <-
                                max_lag
                            data_set_results$filtering[results_row_ind] <-
                                filtering
                            data_set_results$seasonality[results_row_ind] <-
                                seasonality
                            data_set_results$bw_parameterization[results_row_ind] <-
                                bw_parameterization
                            data_set_results$prediction_time[results_row_ind] <-
                                data$time[prediction_time_ind]
                            
                            ## Load kcde_fit object.  Estimation was performed previously.
                            case_descriptor <- paste0(
                                data_set,
                                "-prediction_horizon_", prediction_horizon,
                                "-max_lag_", max_lag,
                                "-filtering_", filtering,
                                "-seasonality_", seasonality,
                                "-bw_parameterization_", bw_parameterization
                            )
                            
                            kcde_fit <- readRDS(file.path(results_path,
                                    paste0("kcde_fit-", case_descriptor, ".rds")))
                            
                            ## Get index of analysis time in data set
                            ## (time from which we predict forward)
                            analysis_time_ind <- prediction_time_ind - prediction_horizon
                            
                            ## Compute log score
                            observed_prediction_target <-
                                data[prediction_time_ind, prediction_target_var, drop = FALSE]
                            colnames(observed_prediction_target) <-
                                paste0(prediction_target_var, "_horizon", prediction_horizon)
                            data_set_results$log_score[results_row_ind] <-
                                kcde_predict(
                                    kcde_fit = kcde_fit,
                                    prediction_data =
                                        data[seq_len(analysis_time_ind), , drop = FALSE],
                                    leading_rows_to_drop = 0L,
                                    trailing_rows_to_drop = 0L,
                                    additional_training_rows_to_drop = NULL,
                                    prediction_type = "distribution",
                                    prediction_test_lead_obs = observed_prediction_target,
                                    log = TRUE
                                )
                            
                            ## Compute point prediction and interval predictions -- quantiles
                            data_set_results[results_row_ind,
                                c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                                kcde_predict(
                                    p = c(0.5, 0.025, 0.975, 0.25, 0.75),
                                    n = 100000,
                                    kcde_fit = kcde_fit,
                                    prediction_data =
                                        data[seq_len(analysis_time_ind), , drop = FALSE],
                                    leading_rows_to_drop = 0L,
                                    trailing_rows_to_drop = 0L,
                                    additional_training_rows_to_drop = NULL,
                                    prediction_type = "quantile",
                                    log = TRUE
                                )
                            
                            ## Correction by subtracting 1 for Dengue data set
                            if(identical(data_set, "dengue_sj")) {
                                data_set_results[results_row_ind,
                                    c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                                    data_set_results[results_row_ind,
                                        c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] - 1L
                            }
                            
                            ## Compute absolute error of point prediction
                            data_set_results$AE[results_row_ind] <-
                                abs(data_set_results$pt_pred[results_row_ind] -
                                    observed_prediction_target)
                            
                            ## Increment results row
                            results_row_ind <- results_row_ind + 1L
                        } # prediction_time_ind
                    } # bw_parameterization
                } # seasonality
            } # filtering
        } # max_lag
    } # prediction_horizon
    
    ## Save results for the given data set
    saveRDS(data_set_results, file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/results",
        data_set,
        "prediction-results/kcde-predictions.rds"))
} # data_set
