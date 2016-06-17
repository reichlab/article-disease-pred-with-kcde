library(kcde)
library(pdtmvn)
library(plyr)
library(dplyr)
library(lubridate)
library(doMC)

options(error = recover)

all_data_sets <- c("ili_national", "dengue_sj")
all_prediction_horizons <- seq_len(52) # make predictions at every horizon from 1 to 52 weeks ahead
all_max_lags <- 1L # use incidence at times t^* and t^* - 1 to predict incidence after t^*
all_max_seasonal_lags <- 0L # not used
all_filtering_values <- FALSE # not used
all_differencing_values <- FALSE # not used
all_seasonality_values <- c(FALSE, TRUE) # specifications without and with periodic kernel
all_bw_parameterizations <- c("diagonal", "full") # specifications with diagonal and full bandwidth

## Test values for debugging
#all_data_sets <- c("ili_national")
#all_prediction_horizons <- 4L
#all_max_lags <- 1L
#all_filtering_values <- c(FALSE)
#all_differencing_values <- c(FALSE)
#all_seasonality_values <- c(TRUE)
#all_bw_parameterizations <- c("diagonal")

num_cores <- 3L
registerDoMC(cores = num_cores)

for(data_set in all_data_sets) {
    ## Set path where fit object is stored
    results_path <- file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results")
    
    if(identical(data_set, "ili_national")) {
        ## Load data for nationally reported influenza like illness
        usflu <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/usflu.csv")
        
#            ## This is how I originally got the data -- have saved it to
#            ## csv for the purposes of stable access going forward.
#            library(cdcfluview)
#            usflu <- get_flu_data("national", "ilinet", years=1997:2014)
        
        ## A little cleanup
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
        data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
        
        ## Season column: for example, weeks of 2010 up through and including week 30 get season 2009/2010;
        ## weeks after week 30 get season 2010/2011
        data$season <- ifelse(
            data$week <= 30,
            paste0(data$year - 1, "/", data$year),
            paste0(data$year, "/", data$year + 1)
        )
        
        ## Season week column: week number within season
        data$season_week <- sapply(seq_len(nrow(data)), function(row_ind) {
                sum(data$season == data$season[row_ind] & data$time_index <= data$time_index[row_ind])
            })
    } else if(identical(data_set, "dengue_sj")) {
        ## Load data for Dengue fever in San Juan
        data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Testing_Data.csv")
        
        ## Form variable with total cases + 0.5 which can be logged, and its seasonally lagged ratio
        data$total_cases_plus_0.5 <- data$total_cases + 0.5
        
        ## Row indices in data corresponding to times at which we want to make a prediction
        prediction_time_inds <- which(data$season %in% paste0(2009:2012, "/", 2010:2013))
        
        ## convert dates
        data$time <- ymd(data$week_start_date)
        
        ## Add time_index column.  This is used for calculating the periodic kernel.
        ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
        ## The origin is arbitrary.
        data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
    }
    
    
    data_set_results <- foreach(prediction_horizon=all_prediction_horizons,
        .packages=c("kcde", "pdtmvn", "plyr", "dplyr", "lubridate"),
        .combine="rbind") %dopar% {
        results_row_ind <- 1L
        
        ## Allocate data frame to store results for this prediction horizon
        num_rows <- length(all_max_lags) *
            length(all_max_seasonal_lags) *
            length(all_filtering_values) *
            length(all_differencing_values) *
            length(all_seasonality_values) *
            length(all_bw_parameterizations) *
            length(prediction_time_inds)
        
        ph_results <- data.frame(data_set = data_set,
            prediction_horizon = rep(NA_integer_, num_rows),
            max_lag = rep(NA_integer_, num_rows),
            max_seasonal_lag = rep(NA_integer_, num_rows),
            filtering = rep(NA, num_rows),
            differencing = rep(NA, num_rows),
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
        class(ph_results$prediction_time) <- class(data$time)
        
        for(max_lag in all_max_lags) {
            for(max_seasonal_lag in all_max_seasonal_lags) {
                for(filtering in all_filtering_values) {
                    for(differencing in all_differencing_values) {
                        ## Set prediction target var
                        if(identical(data_set, "ili_national")) {
                            if(differencing) {
                                prediction_target_var <- "weighted_ili_ratio"
                                orig_prediction_target_var <- "weighted_ili"
                            } else {
                                prediction_target_var <- "weighted_ili"
                            }
                        } else if(identical(data_set, "dengue_sj")) {
                            if(differencing) {
                                prediction_target_var <- "total_cases_plus_0.5_ratio"
                                orig_prediction_target_var <- "total_cases_plus_0.5"
                            } else {
                                prediction_target_var <- "total_cases_plus_0.5"
                            }
                        }
                        
                        for(seasonality in all_seasonality_values) {
                            for(bw_parameterization in all_bw_parameterizations) {
                                for(prediction_time_ind in prediction_time_inds) {
                                    ## Set values describing case in ph_results
                                    ph_results$prediction_horizon[results_row_ind] <-
                                        prediction_horizon
                                    ph_results$max_lag[results_row_ind] <-
                                        max_lag
                                    ph_results$max_seasonal_lag[results_row_ind] <-
                                        max_seasonal_lag
                                    ph_results$filtering[results_row_ind] <-
                                        filtering
                                    ph_results$differencing[results_row_ind] <-
                                        differencing
                                    ph_results$seasonality[results_row_ind] <-
                                        seasonality
                                    ph_results$bw_parameterization[results_row_ind] <-
                                        bw_parameterization
                                    ph_results$prediction_time[results_row_ind] <-
                                        data$time[prediction_time_ind]
                                    
                                    ## Load kcde_fit object.  Estimation was performed previously.
                                    case_descriptor <- paste0(
                                        data_set,
                                        "-prediction_horizon_", prediction_horizon,
                                        "-max_lag_", max_lag,
                                        "-max_seasonal_lag_", max_seasonal_lag,
                                        "-filtering_", filtering,
                                        "-differencing_", differencing,
                                        "-seasonality_", seasonality,
                                        "-bw_parameterization_", bw_parameterization
                                    )
                                    
                                    kcde_fit_file_path <- file.path(results_path,
                                        paste0("kcde_fit-", case_descriptor, ".rds"))
                                    kcde_fit <- readRDS(kcde_fit_file_path)
                                    
                                    
                                    ## If Dengue fit, fix a bug in the in_range function for discrete variables 
                                    ## that was supplied at time of estimation.  This function was not used in
                                    ## estimation, but is required here.  Original definition referenced a function
                                    ## that may not be available now.
                                    if(identical(data_set, "dengue_sj")) {
                                        for(kernel_component_ind in seq_along(kcde_fit$kcde_control$kernel_components)) {
                                            kernel_component <- kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]
                                            
                                            if(!is.null(kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$theta_fixed$discrete_var_range_fns)) {
                                                for(var_ind in seq_along(kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$theta_fixed$discrete_var_range_fns)) {
                                                    kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$theta_fixed$discrete_var_range_fns[[var_ind]]$in_range <-
                                                        function(x, tolerance = .Machine$double.eps^0.5) {
                                                        return(sapply(x,
                                                                function(x_i) {
                                                                    return(isTRUE(all.equal(x_i, 
                                                                                kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$theta_fixed$discrete_var_range_fns[[var_ind]]$discretizer(x_i))))
                                                                }))
                                                    }
                                                }
                                                
                                                kcde_fit$theta_hat[[kernel_component_ind]]$discrete_var_range_fns <-
                                                    kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$theta_fixed$discrete_var_range_fns
                                            }
                                        }
                                    }
                                    
                                    ## fix rkernel_fn for pdtmvn-based kernel functions
                                    ## I supplied a buggy version of this in the call to the estimation routine
                                    ## that did not ensure that variables were supplied in a consistent order.
                                    ## This did not affect estimation as rkernel_fn is not called there
                                    ## But it does need to be fixed here.
                                    for(kernel_component_ind in seq(from = (as.logical(seasonality) + 1), to = length(kcde_fit$kcde_control$kernel_components))) {
                                        kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$rkernel_fn <-
                                            function(n,
                                                conditioning_obs,
                                                center,
                                                bw,
                                                bw_continuous,
                                                conditional_bw_discrete,
                                                conditional_center_discrete_offset_multiplier,
                                                continuous_vars,
                                                discrete_vars,
                                                continuous_var_col_inds,
                                                discrete_var_col_inds,
                                                discrete_var_range_fns,
                                                lower,
                                                upper,
                                                x_names,
                                                ...) {
                                            if(missing(conditioning_obs) || is.null(conditioning_obs)) {
                                                log_conditioning_obs <- NULL
                                            } else {
                                                log_conditioning_obs <- log(conditioning_obs)
                                            }
                                            if(missing(bw_continuous)) {
                                                bw_continuous <- NULL
                                            }
                                            if(missing(conditional_bw_discrete)) {
                                                conditional_bw_discrete <- NULL
                                            }
                                            if(missing(conditional_center_discrete_offset_multiplier)) {
                                                conditional_center_discrete_offset_multiplier <- NULL
                                            }
                                            
                                            ## center parameter of pdtmvn_kernel is mean of log
                                            ## mode of resulting log-normal distribution is
                                            ## mode = exp(mu - bw %*% 1) (where 1 is a column vector of 1s)
                                            ## therefore mu = log(mode) + bw %*% 1
                                            reduced_x_names <- names(center)
                                            inds_x_vars_in_orig_vars <- which(x_names %in% reduced_x_names)
                                            x_names_for_call <- x_names[inds_x_vars_in_orig_vars]
                                            
                                            mean_offset <- apply(bw, 1, sum)[x_names %in% colnames(center)]
                                            
                                            return(exp(rpdtmvn_kernel(n = n,
                                                        conditioning_obs = log_conditioning_obs,
                                                        center = sweep(log(center)[, x_names_for_call, drop = FALSE], 2, mean_offset, `+`),
                                                        bw = bw,
                                                        bw_continuous = bw_continuous,
                                                        conditional_bw_discrete = conditional_bw_discrete,
                                                        conditional_center_discrete_offset_multiplier = conditional_center_discrete_offset_multiplier,
                                                        continuous_vars = continuous_vars,
                                                        discrete_vars = discrete_vars,
                                                        continuous_var_col_inds = continuous_var_col_inds,
                                                        discrete_var_col_inds = discrete_var_col_inds,
                                                        discrete_var_range_fns = discrete_var_range_fns,
                                                        lower = lower,
                                                        upper = upper,
                                                        x_names = x_names)[, reduced_x_names, drop = FALSE]))
                                        }
                                    }
                                    
                                    
                                    ## Get index of analysis time in data set
                                    ## (time from which we predict forward)
                                    analysis_time_ind <- prediction_time_ind - prediction_horizon
                                    
                                    ## Compute log score
                                    observed_prediction_target <-
                                        data[prediction_time_ind, prediction_target_var, drop = FALSE]
                                    colnames(observed_prediction_target) <-
                                        paste0(prediction_target_var, "_horizon", prediction_horizon)
                                    ph_results$log_score[results_row_ind] <-
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
                                    
                                    if(differencing) {
                                        ph_results$log_score[results_row_ind] <-
                                            ph_results$log_score[results_row_ind] -
                                            (abs(data[analysis_time_ind - 52, orig_prediction_target_var]))
                                    }
                                    
                                    ## Compute point prediction and interval predictions -- quantiles
                                    ph_results[results_row_ind,
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
                                    
                                    if(differencing) {
                                        ph_results[results_row_ind,
                                            c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                                            ph_results[results_row_ind,
                                                c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] * 
                                            data[prediction_time_ind - 52, prediction_target_var]
                                    }
                                    
                                    ## Correction by subtracting 0.5 for Dengue data set
                                    if(identical(data_set, "dengue_sj")) {
                                        ph_results[results_row_ind,
                                            c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                                            ph_results[results_row_ind,
                                                c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] - 0.5
                                    }
                                    
                                    ## Compute absolute error of point prediction
                                    ph_results$AE[results_row_ind] <-
                                        abs(ph_results$pt_pred[results_row_ind] -
                                                observed_prediction_target)
                                    
                                    ## Increment results row
                                    results_row_ind <- results_row_ind + 1L
                                } # prediction_time_ind
                            } # bw_parameterization
                        } # seasonality
                    } # differencing
                } # filtering
            } # max_seasonal_lag
        } # max_lag
        
        saveRDS(ph_results, file = file.path(
                "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
                data_set,
                "prediction-results",
                paste0("kcde-predictions-ph_", prediction_horizon, ".rds")))
        
        return(ph_results)
    } # prediction_horizon -- in foreach/dopar statement
    
    ## Save results for the given data set
    saveRDS(data_set_results, file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "prediction-results/kcde-predictions.rds"))
} # data_set
