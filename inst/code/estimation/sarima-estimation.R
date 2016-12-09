library(lubridate)
library(forecast)

all_data_sets <- c("ili_national", "dengue_sj")
all_data_sets <- "dengue_sj"

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
        
        ## Subset data to do estimation using only data up through 2010
        ## 2011 - 2014 are held out for evaluating performance.
        data <- data[data$year <= 2010, , drop = FALSE]
        
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
        data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Training_Data.csv")
        
        ## Restrict to data from 1990/1991 through 2008/2009 seasons
        train_seasons <- paste0(1990:2008, "/", 1991:2009)
        data <- data[data$season %in% train_seasons, ]
        
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
    
    seasonally_differenced_log_sarima_fit <-
        auto.arima(seasonally_differenced_log_prediction_target)
    
    saveRDS(seasonally_differenced_log_sarima_fit,
        file = file.path(
            "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
            data_set,
            "estimation-results/sarima-fit.rds"))
}


predictions_df <- data.frame(ph=rep(seq_len(52), times = 5 * 52),
	last_obs_season=rep(c("2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013"), each = 52, times = 52),
	last_obs_week=rep(seq_len(52) - 1, each = 52 * 5),
	model="sarima",
	stringsAsFactors=FALSE)
predictions_df$prediction_season <- predictions_df$last_obs_season
predictions_df$prediction_week <- predictions_df$last_obs_week + predictions_df$ph

inds_last_obs_season_prev_year <- which(predictions_df$last_obs_week == 0)
predictions_df$last_obs_season[inds_last_obs_season_prev_year] <- 
	sapply(predictions_df$last_obs_season[inds_last_obs_season_prev_year],
		function(next_season) {
			start_year <- as.integer(substr(next_season, 1, 4)) - 1L
			paste0(start_year, "/", start_year + 1)
		}
	)
predictions_df$last_obs_week[inds_last_obs_season_prev_year] <- 52L

inds_prediction_season_next_year <- which(predictions_df$prediction_week > 52)
predictions_df$prediction_season[inds_prediction_season_next_year] <- 
	sapply(predictions_df$prediction_season[inds_prediction_season_next_year],
		function(next_season) {
			start_year <- as.integer(substr(next_season, 1, 4)) + 1L
			paste0(start_year, "/", start_year + 1)
		}
	)
predictions_df$prediction_week[inds_prediction_season_next_year] <-
	predictions_df$prediction_week[inds_prediction_season_next_year] - 52L

predictions_df$log_score <- NA
predictions_df$prediction <- NA
predictions_df$AE <- NA
predictions_df$predictive_50pct_lb <- NA
predictions_df$predictive_50pct_ub <- NA
predictions_df$predictive_90pct_lb <- NA
predictions_df$predictive_90pct_ub <- NA
predictions_df$week_start_date <- San_Juan_test$week_start_date[1]



#sj_log_sarima_fit <- sj_log_sarima_fit_manual1  ## very bad
#sj_log_sarima_fit <- sj_log_sarima_fit_manual2  ## quite bad
#sj_log_sarima_fit <- sj_log_sarima_fit_manual3  ## very bad

sarima_inds <- which(predictions_df$model == "sarima")
#sarima_inds <- which(predictions_df$model == "sarima" & predictions_df$ph %in% c(1, 13, 26, 39, 52))
for(predictions_df_row_ind in sarima_inds) {
	ph <- as.numeric(predictions_df$ph[predictions_df_row_ind])
	last_obs_ind <- which(San_Juan_test$season == predictions_df$last_obs_season[predictions_df_row_ind] &
		San_Juan_test$season_week == predictions_df$last_obs_week[predictions_df_row_ind])
	
	predictions_df$week_start_date[predictions_df_row_ind] <- San_Juan_test$week_start_date[last_obs_ind + as.numeric(ph)]
	
	
#	new_data <- ts(log(San_Juan_test$total_cases[seq_len(last_obs_ind)] + 1), frequency = 52)
#	updated_sj_log_sarima_fit <- Arima(new_data, model = sj_log_sarima_fit)


	new_data <- log(San_Juan_test$total_cases[seq_len(last_obs_ind)] + 1)
	seasonal_diff_new_data <- ts(new_data[seq(from = 53, to = length(new_data))] -
	    new_data[seq(from = 1, to = length(new_data) - 52)], frequency = 52)
	updated_sj_log_sarima_fit <- Arima(seasonal_diff_new_data, model = sj_seasonally_differenced_log_sarima_fit)
	
	predict_result <- predict(updated_sj_log_sarima_fit, n.ahead = ph)
	
#	predictive_log_mean <- as.numeric(predict_result$pred[ph])
#	predictions_df$prediction[predictions_df_row_ind] <- exp(predictive_log_mean) - 1
	
	predictive_log_mean <- as.numeric(predict_result$pred[ph]) + new_data[last_obs_ind + ph - 52]
	predictions_df$prediction[predictions_df_row_ind] <- exp(as.numeric(predict_result$pred[ph]) + new_data[last_obs_ind + ph - 52]) - 1


	predictions_df$AE[predictions_df_row_ind] <- abs(predictions_df$prediction[predictions_df_row_ind] - San_Juan_test$total_cases[last_obs_ind + ph])
	predictions_df$log_score[predictions_df_row_ind] <- dlnorm(San_Juan_test$total_cases[last_obs_ind + ph] + 1,
		meanlog = predictive_log_mean,
		sdlog = as.numeric(predict_result$se[ph]),
		log = TRUE)
	temp <- qlnorm(c(0.05, 0.25, 0.75, 0.95),
		meanlog = predictive_log_mean,
		sdlog = as.numeric(predict_result$se[ph]))
	predictions_df[predictions_df_row_ind, c("predictive_90pct_lb", "predictive_50pct_lb", "predictive_50pct_ub", "predictive_90pct_ub")] <-
		temp - 1
}

predictions_df <- predictions_df[predictions_df$week_start_date %in% San_Juan_test$week_start_date[San_Juan_test$season %in% c("2009/2010", "2010/2011", "2011/2012", "2012/2013")], ]

predictions_df$ph <- as.factor(predictions_df$ph)

sarima_predictions_df <- predictions_df
save(sarima_predictions_df, file = "/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/sarima-predictions.Rdata")

library(ggplot2)
ggplot() +
	geom_line(aes(x = week_start_date, y = total_cases), data = San_Juan_test) +
	geom_line(aes(x = week_start_date, y = prediction, colour = ph), data = predictions_df[predictions_df$ph %in% c(1, 13, 26, 39, 52), , drop = FALSE]) +
	theme_bw()


