library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(surveillance)
library(mvtnorm)
library(doMC)

set.seed(760933)

data_set <- "dengue_sj"

n_sims <- 10000

prediction_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
    data_set,
    "prediction-results")



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

seasonally_differenced_log_prediction_target <-
    ts(log_prediction_target[seq(from = 53, to = length(log_prediction_target))] -
            log_prediction_target[seq(from = 1, to = length(log_prediction_target) - 52)],
        frequency = 52)

## load surveillance fits and choose best one
surveillance_fits <- readRDS(file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results/surveillance-fits.rds"))

best_spec_ind <- which.min(surveillance_fits$model_specifications$mean_log_score)
surveillance_fit <- surveillance_fits$model_fits[[best_spec_ind]]



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
    
    observed_peak_height <- which(
        dengue_incidence_bins$lower <= observed_peak_height &
            dengue_incidence_bins$upper > observed_peak_height)
    
    for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
        analysis_time_ind <- which(data$season == analysis_time_season &
                data$season_week == analysis_time_season_week)
        
        ## simulate trajectories
        max_prediction_horizon <-
          last_analysis_time_season_week + 1 -
          analysis_time_season_week
        
        trajectory_samples <- t(simulate(surveillance_fit,
            nsim = n_sims,
            y.start = data$total_cases[analysis_time_ind],
            subset = seq(from = analysis_time_ind + 1, to = analysis_time_ind + max_prediction_horizon),
            simplify = TRUE)[, 1, ])

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
        paste0("peak-week-surveillance-", data_set, ".rds")))
