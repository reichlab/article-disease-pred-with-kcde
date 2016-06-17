library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(forecast)
library(mvtnorm)
library(doMC)


data_set <- "ili_national"
data_set <- "dengue_sj"

n_sims <- 10000

prediction_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
    data_set,
    "prediction-results")

#sample_predictive_trajectories_arima <- function(n, n.ahead, arima_fit) {
#    results <- matrix(NA, nrow = n, ncol = n.ahead)
#    for(i in seq_len(n)) {
#        updated_arima_fit <- arima_fit
#        for(j in seq_len(n.ahead)) {
#            predict_result <- predict(updated_arima_fit, n.ahead = 1)
#            results[i, j] <- rnorm(1, mean = predict_result$pred[1], sd = predict_result$se[1])
#            updated_arima_fit$x <-
#                c(updated_arima_fit$x, results[i, j])
#        }
#    }
#    
#    return(results)
#}

## This is directly taken from the forecast.Arima function from the forecast package,
## but I've forced bootstrap = TRUE and return the simulated trajectories which were
## not returned in the original function definition.
sample_predictive_trajectories_arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 
        10), level = c(80, 95), fan = FALSE, xreg = NULL, lambda = object$lambda, 
        npaths = 5000, ...) 
{
    sim <- matrix(NA, nrow = npaths, ncol = h)
    
    for (i in 1:npaths) {
        sim[i, ] <- simulate.Arima(object,
            nsim = h)
    }
    
    return(sim)
}




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
    
    
    season_length <- 33L
    analysis_seasons <- c("2011/2012", "2012/2013", "2013/2014")
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    
    prediction_target_var <- "weighted_ili"
    
    log_prediction_target <- log(data[, prediction_target_var])
    
    ili_incidence_bins <- data.frame(
        lower = seq(from = 0, to = 13, by = 0.5),
        upper = c(seq(from = 0.5, to = 13, by = 0.5), Inf))
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
    
    season_length <- 52L
    analysis_seasons <- paste0(2009:2012, "/", 2010:2013)
    first_analysis_time_season_week <- 1 # == week 40 of year
    last_analysis_time_season_week <- 51 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    
    prediction_target_var <- "total_cases_plus_1"
    
    log_prediction_target <- log(data[, prediction_target_var])
    
    dengue_incidence_bins <- data.frame(
        lower = seq(from = 0, to = 500, by = 50),
        upper = c(seq(from = 50, to = 500, by = 50), Inf))
}

seasonally_differenced_log_prediction_target <-
    ts(log_prediction_target[seq(from = 53, to = length(log_prediction_target))] -
            log_prediction_target[seq(from = 1, to = length(log_prediction_target) - 52)],
        frequency = 52)

seasonally_differenced_log_sarima_fit <- readRDS(file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results/sarima-fit.rds"))

results <- cbind(
        expand.grid(analysis_seasons,
            seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1),
            stringsAsFactors = FALSE),
        matrix(NA,
            nrow = length(analysis_seasons) * (last_analysis_time_season_week - first_analysis_time_season_week),
            ncol = 3 * n_sims + 2)
    ) %>%
    `colnames<-`(c("analysis_time_season",
            "analysis_time_season_week",
            "peak_week_log_score",
            "peak_height_log_score",
            paste0("peak_week_", seq_len(n_sims)),
            paste0("peak_height_", seq_len(n_sims)),
            paste0("unbinned_peak_height_", seq_len(n_sims))))

## generate peak week timing and height estimates
for(analysis_time_season in analysis_seasons) {
    ## get observed quantities, for computing log score
    observed_peak_height <- max(data[data$season == analysis_time_season, prediction_target_var])
    observed_peak_week_ind <- which((data$season == analysis_time_season) &
            (data[, prediction_target_var] == observed_peak_height))
    observed_peak_week <- data$season_week[observed_peak_week_ind]
    
    if(identical(data_set, "ili_national")) {
        observed_peak_height <- which(
            ili_incidence_bins$lower <= observed_peak_height &
                ili_incidence_bins$upper > observed_peak_height)
    } else if(identical(data_set, "dengue_sj")) {
        observed_peak_height <- which(
            dengue_incidence_bins$lower <= observed_peak_height &
                dengue_incidence_bins$upper > observed_peak_height)
    }
    
    for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
        analysis_time_ind <- which(data$season == analysis_time_season &
                data$season_week == analysis_time_season_week)
        
        ## Update anova fit object with seasonally differenced data up
        ## through analysis_time_ind
        new_data <- seasonally_differenced_log_prediction_target[
            seq_len(max(0, analysis_time_ind - 52))]
        last_na_ind <- max(c(0, max(which(is.na(new_data)))))
        new_data <- new_data[seq(from = last_na_ind + 1, to = length(new_data))]
        updated_log_sarima_fit <- Arima(
            new_data,
            model = seasonally_differenced_log_sarima_fit)
        
        ## simulate n_sims trajectories recursively from sarima
        max_prediction_horizon <-
            last_analysis_time_season_week + 1 -
            analysis_time_season_week
        trajectory_samples <- sample_predictive_trajectories_arima(
            updated_log_sarima_fit,
            h = max_prediction_horizon,
            npaths = n_sims)
        
        ## Sampled trajectories are of seasonally differenced log incidence
        ## Get to trajectories for originally observed incidence by
        ## adding seasonal lag of incidence and exponentiating
        for(prediction_horizon in seq_len(max_prediction_horizon)) {
            trajectory_samples[, prediction_horizon] <- trajectory_samples[, prediction_horizon] +
                log_prediction_target[analysis_time_ind + prediction_horizon - 52]
        }
        trajectory_samples <- exp(trajectory_samples)
        
        ## Augment trajectory samples with previous observed incidence values
        season_start_ind <- which(data$season == analysis_time_season &
                data$season_week == 1)
        if(season_start_ind < analysis_time_ind) {
            trajectory_samples <- cbind(
                matrix(
                    rep(data[seq(from = season_start_ind, to = analysis_time_ind), prediction_target_var], each = n_sims),
                    nrow = n_sims
                ),
                trajectory_samples
            )
        }
        
        ## Get peak week and height at peak week for each simulated trajectory 
        results_save_row <- which(results$analysis_time_season == analysis_time_season &
                results$analysis_time_season_week == analysis_time_season_week)
        
        peak_week_by_sim_ind <- apply(trajectory_samples, 1, which.max)
        results[results_save_row, paste0("peak_week_", seq_len(n_sims))] <-
            peak_week_by_sim_ind
        
        peak_week_height_by_sim_ind <- trajectory_samples[cbind(seq_len(n_sims), peak_week_by_sim_ind)]
        results[results_save_row, paste0("unbinned_peak_height_", seq_len(n_sims))] <-
            peak_week_height_by_sim_ind
        
        if(identical(data_set, "ili_national")) {
            peak_week_height_by_sim_ind <- sapply(peak_week_height_by_sim_ind,
                function(height) {
                    which(ili_incidence_bins$lower <= height &
                            ili_incidence_bins$upper > height)
                })
        } else if(identical(data_set, "dengue_sj")) {
            peak_week_height_by_sim_ind <- sapply(peak_week_height_by_sim_ind,
                function(height) {
                    which(dengue_incidence_bins$lower <= height &
                            dengue_incidence_bins$upper > height)
                })
        }
        results[results_save_row, paste0("peak_height_", seq_len(n_sims))] <-
            peak_week_height_by_sim_ind
        
        ## Get log scores
        results[results_save_row, "peak_week_log_score"] <- log(sum(peak_week_by_sim_ind == observed_peak_week)) - log(n_sims)
        results[results_save_row, "peak_height_log_score"] <- log(sum(peak_week_height_by_sim_ind == observed_peak_height)) - log(n_sims)
    }
}

saveRDS(results,
    file.path(prediction_save_path,
        paste0("peak-week-sarima-", data_set, ".rds")))


## Visual check
#library(plyr)
#library(dplyr)
#library(tidyr)
#library(ggplot2)
#
#prev_preds <- readRDS("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/ili_national/prediction-results/sarima-predictions.rds")
#inds_keep <- sapply(1:33, function(ph) {
#    which(prev_preds$prediction_time == data$time[analysis_time_ind + ph] & prev_preds$prediction_horizon == ph)
#})
#
#ts_for_plot <- trajectory_samples %>%
#    as.data.frame() %>%
#    `colnames<-`(as.character(seq(from = -1 * (analysis_time_season_week - 1), length = ncol(trajectory_samples)))) %>%
#    mutate(sim_ind = seq_len(nrow(trajectory_samples))) %>%
#    gather_("prediction_horizon", "sim_incidence", as.character(seq(from = -1 * (analysis_time_season_week - 1), length = ncol(trajectory_samples)))) %>%
#    mutate(prediction_time = data$time[analysis_time_ind + as.integer(prediction_horizon)])
#
#ggplot() +
#    geom_line(aes(x = as.Date(time), y = weighted_ili), data = data[seq(from = analysis_time_ind - 10, to = analysis_time_ind + 52), ]) +
#    geom_line(aes(x = as.Date(prediction_time), y = sim_incidence, group = sim_ind), colour = "blue", alpha = 0.5, data = ts_for_plot) +
#    geom_line(aes(x = as.Date(prediction_time), y = pt_pred), colour = "red", data = prev_preds[inds_keep, ]) +
##    geom_hline(yintercept = 4.25906) +
#    geom_point(aes(x = as.Date(data$time[analysis_time_ind + peak_week_by_sim_ind - (analysis_time_season_week)]), y = peak_week_height_by_sim_ind)) +
#    scale_y_log10() +
#    theme_bw()
