library(forecast)
library(cdcfluview)
library(plyr)
library(dplyr)
library(kcde)
library(lubridate)

all_data_sets <- c("ili_national", "dengue_sj")
all_prediction_horizons <- 1:52
all_prediction_statistics <- c("log_score",
    "pt_pred",
    "AE",
    "interval_pred_lb_95",
    "interval_pred_ub_95",
    "interval_pred_lb_50",
    "interval_pred_ub_50")

for(data_set in all_data_sets) {
    ### Load data set and set variables describing how the fit is performed
    if(identical(data_set, "ili_national")) {
        ## Load data for nationally reported influenza like illness
        usflu <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/usflu.csv")
        
#            ## This is how I originally got the data -- have saved it to
#            ## csv for the purposes of stable access going forward.
#            library(cdcfluview)
#            usflu <- get_flu_data("national", "ilinet", years=1997:2014)
        
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
        
        prediction_target_var <- "weighted_ili"
        
        log_prediction_target <- log(data[, prediction_target_var])
    } else if(identical(data_set, "dengue_sj")) {
        ## Load data for Dengue fever in San Juan
        data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Testing_Data.csv")
        
        ## Row indices in data corresponding to times at which we want to make a prediction
        prediction_time_inds <- which(data$season %in% paste0(2009:2012, "/", 2010:2013))
        
        ## Form variable with total cases + 1 which can be logged
        data$total_cases_plus_1 <- data$total_cases + 1
        
        ## convert dates
        data$time <- ymd(data$week_start_date)
        
        ## Add time_index column.  This is used for calculating the periodic kernel.
        ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
        ## The origin is arbitrary.
        data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
        
        prediction_target_var <- "total_cases_plus_1"
        
        log_prediction_target <- log(data[, prediction_target_var])
    }
    
    seasonally_differenced_log_prediction_target <-
        ts(log_prediction_target[seq(from = 53, to = length(log_prediction_target))] -
                log_prediction_target[seq(from = 1, to = length(log_prediction_target) - 52)],
            frequency = 52)
    
    seasonally_differenced_log_sarima_fit <- readRDS(file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results/sarima-fit.rds"))
    
    
    num_rows <- length(all_prediction_horizons) *
        length(prediction_time_inds)
    
    data_set_results <- data.frame(data_set = data_set,
        model = "SARIMA",
        prediction_horizon = rep(NA_integer_, num_rows),
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
        for(prediction_time_ind in prediction_time_inds) {
            ## Set values describing case in data_set_results
            data_set_results$prediction_horizon[results_row_ind] <-
                prediction_horizon
            data_set_results$prediction_time[results_row_ind] <-
                data$time[prediction_time_ind]
            
            ## Get index of analysis time in data set
            ## (time from which we predict forward)
            analysis_time_ind <- prediction_time_ind - prediction_horizon
            
            ## Observed value at prediction time -- used in calculating log
            ## score and absolute error
            observed_prediction_target <-
                data[prediction_time_ind, prediction_target_var]
            
            ## Update anova fit object with seasonally differenced data up
            ## through analysis_time_ind
            new_data <- seasonally_differenced_log_prediction_target[
                seq_len(max(0, analysis_time_ind - 52))]
            updated_log_sarima_fit <- Arima(
                new_data,
                model = seasonally_differenced_log_sarima_fit)
            
            ## Get mean and variance of predictive distribution on log scale
            ## Call to predict gets an object with mean mu* and sd sd* parameters
            ## for the (normally distributed) predictive distribution for
            ## log(incidence_{t* + prediction_horizon} / incidence_{t* + prediction_horizon - 52})
            ## This means that the predictive distribution for incidence_{t* + prediction_horizon} is
            ## lognorm(mu* + log(incidence_{t* + prediction_horizon - 52}), sd*)
            predict_result <- predict(updated_log_sarima_fit, n.ahead = prediction_horizon)
            predictive_log_mean <- as.numeric(predict_result$pred[prediction_horizon]) +
                log_prediction_target[analysis_time_ind + prediction_horizon - 52]
            predictive_log_sd <- as.numeric(predict_result$se[prediction_horizon])
            
            ## Compute log score of distribution prediction
            if(identical(data_set, "ili_national")) {
                data_set_results$log_score[results_row_ind] <-
                    dlnorm(observed_prediction_target,
                        meanlog = predictive_log_mean,
                        sdlog = predictive_log_sd,
                        log = TRUE)
            } else if(identical(data_set, "dengue_sj")) {
                data_set_results$log_score[results_row_ind] <-
                    logspace_sub(plnorm(observed_prediction_target + 0.5,
                        meanlog = predictive_log_mean,
                        sdlog = predictive_log_sd,
                        log = TRUE),
                    plnorm(observed_prediction_target - 0.5,
                        meanlog = predictive_log_mean,
                        sdlog = predictive_log_sd,
                        log = TRUE))
            }
            
            ## Compute point prediction
            data_set_results$pt_pred[results_row_ind] <-
                exp(predictive_log_mean)
            
            ## Compute absolute error of point prediction
            data_set_results$AE[results_row_ind] <-
                abs(data_set_results$pt_pred[results_row_ind] -
                        observed_prediction_target)
            
            ## Correct point prediction by subtracting 1 for Dengue data set
            ## We didn't do this before computing absolute error because
            ## observed_prediction_target has 1 added to it.
            if(identical(data_set, "dengue_sj")) {
                data_set_results$pt_pred[results_row_ind] <-
                    data_set_results$pt_pred[results_row_ind] - 1L
            }
            
            
            ## Compute prediction interval bounds
            data_set_results[results_row_ind,
                c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                    qlnorm(c(0.025, 0.975, 0.25, 0.75),
                        meanlog = predictive_log_mean,
                        sdlog = predictive_log_sd)
            
            ## Correction by subtracting 1 for Dengue data set
            if(identical(data_set, "dengue_sj")) {
                data_set_results[results_row_ind,
                    c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                    data_set_results[results_row_ind,
                        c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] - 1L
            }
            
            ## Increment results row
            results_row_ind <- results_row_ind + 1L
        } # prediction_time_ind
    } # prediction_horizon
    
    ## Save results for the given data set
    saveRDS(data_set_results, file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "prediction-results/sarima-predictions.rds"))
} # data_set
