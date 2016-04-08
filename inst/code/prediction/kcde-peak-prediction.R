library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(copula)
library(mvtnorm)
library(doMC)

## Function borrowed from copula package and modified
makePosDef <- function (mat, delta = 0.001) 
{
    while(min(eigen(mat)$values) < 10^{-6}) {
        decomp <- eigen(mat)
        Lambda <- decomp$values
        Lambda[Lambda < 0] <- delta
        Gamma <- decomp$vectors
        newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
        D <- 1/sqrt(diag(newmat))
        mat <- diag(D) %*% newmat %*% diag(D)
    }
    return(mat)
}


## outer loop runs in parallel over cases defined by combinations of
## data set, max lag, max seasonal lag, filtering, differencing, seasonality, bw_parameterization
data_set <- "ili_national"
if(identical(data_set, "ili_national")) {
    all_max_lags <- as.character(c(1L))
    all_max_seasonal_lags <- as.character(c(0L, 1L))
    all_filtering_values <- c("FALSE")
    all_differencing_values <- c("FALSE", "TRUE")
    all_seasonality_values <- c("FALSE", "TRUE")
    all_bw_parameterizations <- c("diagonal", "full")
    
    incidence_bins <- data.frame(
        lower = seq(from = 0, to = 13, by = 0.5),
        upper = c(seq(from = 0.5, to = 13, by = 0.5), Inf))
    
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
    results_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-results"
    scripts_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-scripts"
} else if(identical(data_set, "dengue_sj")) {
    all_prediction_horizons <- as.character(seq_len(52))
    all_max_lags <- as.character(c(1L))
    all_max_seasonal_lags <- as.character(c(0L, 1L))
    all_filtering_values <- c("FALSE")
#    all_differencing_values <- c("FALSE", "TRUE")
    all_differencing_values <- "FALSE"
    all_seasonality_values <- c("FALSE", "TRUE")
    all_bw_parameterizations <- c("diagonal", "full")
    all_sim_n <- "NA"
    all_sim_families <- "NA"
    all_sim_run_inds <- 1L
    
    incidence_bins <- data.frame(
        lower = seq(from = 0, to = 500, by = 50),
        upper = c(seq(from = 50, to = 500, by = 50), Inf))
    
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
    results_path <- "/home/er71a/kcde-applied-paper/R/application-dengue/estimation-results"
    scripts_path <- "/home/er71a/kcde-applied-paper/R/application-dengue/estimation-scripts"
} else if(identical(data_set, "sim")) {
    all_prediction_horizons <- "0"
    all_max_lags <- "0"
    all_max_seasonal_lags <- "0"
    all_filtering_values <- "FALSE"
    all_differencing_values <- c("FALSE")
    all_seasonality_values <- c("FALSE")
    
    all_bw_parameterizations <- c("diagonal", "full")
#        all_bw_parameterizations <- c("full")
#        all_bw_parameterizations <- c("diagonal")
    all_sim_n <- c("100", "1000")
#        all_sim_n <- c("100")
#        all_sim_families <- c("bivariate-B-discretized",
#            "bivariate-C-discretized",
#        all_sim_families <- c("multivariate-2d-discretized",
#            "multivariate-4d-discretized", 
#            "multivariate-6d-discretized")
    all_sim_families <- c("multivariate-2d-discretized",
        "multivariate-4d-discretized")
#        all_sim_families <- c("multivariate-4d-discretized")
#        all_sim_run_inds <- seq(from = 1, to = 100)
#        all_sim_run_inds <- 9L
    all_sim_run_inds <- c(seq(from = 1, to = 8), seq(from = 10, to = 100))
    
#        all_bw_parameterizations <- c("diagonal")
#        all_sim_n <- c("100")
#        all_sim_families <- c("multivariate-6d-discretized")
#        all_sim_run_inds <- 45L
} else {
    stop("Invalid data set")
}
    

case_definitions <- expand.grid(
        data_set,
        all_max_lags,
        all_max_seasonal_lags,
        all_filtering_values,
        all_differencing_values,
        all_seasonality_values,
        all_bw_parameterizations,
        stringsAsFactors = FALSE) %>%
    `colnames<-`(c("data_set",
            "max_lag",
            "max_seasonal_lag",
            "filtering",
            "differencing",
            "seasonality",
            "bw_parameterization"))

num_cores <- 3L
registerDoMC(cores = num_cores)

junk <- foreach(case_row_ind = seq_len(nrow(case_definitions)),
        .packages = c("kcde", "plyr", "dplyr", "lubridate", "reshape", "copula", "mvtnorm"),
        .combine = "rbind") %dopar% {
    data_set <- case_definitions$data_set[case_row_ind]
    max_lag <- case_definitions$max_lag[case_row_ind]
    max_seasonal_lag <- case_definitions$max_seasonal_lag[case_row_ind]
    filtering <- case_definitions$filtering[case_row_ind]
    differencing <- case_definitions$differencing[case_row_ind]
    seasonality <- case_definitions$seasonality[case_row_ind]
    bw_parameterization <- case_definitions$bw_parameterization[case_row_ind]
    
    n_sims <- 10000
    
    copula_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "copula-estimation-results")
    estimation_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results")
    prediction_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "prediction-results")
    
    case_descriptor <- paste0(
        data_set,
        "-max_lag_", max_lag,
        "-max_seasonal_lag_", max_seasonal_lag,
        "-filtering_", filtering,
        "-differencing_", differencing,
        "-seasonality_", seasonality,
        "-bw_parameterization_", bw_parameterization
    )
    file_name <- paste0("kcde-copula-fits-",
        case_descriptor,
        ".rds")
    copula_fits <- readRDS(file = file.path(copula_save_path, file_name))
    analysis_time_season_week_by_copula_fit <- unlist(lapply(copula_fits,
            function(copula_fit) {copula_fit$analysis_time_season_week}))
    
    
    kcde_fits_by_prediction_horizon <- lapply(seq_len(52),
        function(prediction_horizon) {
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
            readRDS(file.path(estimation_save_path,
                    paste0("kcde_fit-", case_descriptor, ".rds")))
        })
    
    
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
        
        if(differencing) {
            data$weighted_ili_ratio <- data$weighted_ili / lag(data$weighted_ili, 52) 
            orig_prediction_target_var <- "weighted_ili"
            prediction_target_var <- "weighted_ili_ratio"
        } else {
            orig_prediction_target_var <- "weighted_ili"
            prediction_target_var <- "weighted_ili"
        }
        
        analysis_seasons <- c("2011/2012", "2012/2013", "2013/2014")
        first_analysis_time_season_week <- 10 # == week 40 of year
        last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    } else if(identical(data_set, "dengue_sj")) {
        ## Load data for Dengue fever in San Juan
        data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Testing_Data.csv")
        
        ## Form variable with total cases + 1 which can be logged
        data$total_cases_plus_1 <- data$total_cases + 1
        data$total_cases_plus_1_ratio <- data$total_cases_plus_1 / lag(data$total_cases_plus_1, 52)
        
        ## convert dates
        data$time <- ymd(data$week_start_date)
        
        ## Add time_index column.  This is used for calculating the periodic kernel.
        ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
        ## The origin is arbitrary.
        data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
        
        if(differencing) {
            orig_prediction_target_var <- "total_cases_plus_1"
            prediction_target_var <- "total_cases_plus_1_ratio"
        } else {
            orig_prediction_target_var <- "total_cases_plus_1"
            prediction_target_var <- "total_cases_plus_1"
        }
        
        analysis_seasons <- c("2009/2010",
            "2010/2011",
            "2011/2012",
            "2012/2013")
        first_analysis_time_season_week <- 1
        last_analysis_time_season_week <- 51
    }
    
    results <- cbind(
            expand.grid(analysis_seasons,
                seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1),
                stringsAsFactors = FALSE),
            matrix(NA,
                nrow = length(analysis_seasons) * (last_analysis_time_season_week - first_analysis_time_season_week),
                ncol = 2 * n_sims + 2)
        ) %>%
        `colnames<-`(c("analysis_time_season",
                "analysis_time_season_week",
                "peak_week_log_score",
                "peak_height_log_score",
                paste0("peak_week_", seq_len(n_sims)),
                paste0("peak_height_", seq_len(n_sims))))
    
    ## generate peak week timing and height estimates
    for(analysis_time_season in analysis_seasons) {
        ## get observed quantities, for computing log score
        observed_peak_height <- max(data[data$season == analysis_time_season, orig_prediction_target_var])
        observed_peak_week_ind <- which((data$season == analysis_time_season) &
                (data[, orig_prediction_target_var] == observed_peak_height))
        observed_peak_week <- data$season_week[observed_peak_week_ind]
        
        observed_peak_height <- which(
            incidence_bins$lower <= observed_peak_height &
                incidence_bins$upper > observed_peak_height)
        
        for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
            ### simulate from copula that ties the marginal predictive distributions together
            
            ## get the right copula for analysis_time_season_week
            predictive_copula_ind <- which(analysis_time_season_week_by_copula_fit == analysis_time_season_week)
            copula_fit <- copula_fits[[predictive_copula_ind]]$copula_fit
            predictive_copula <- copula_fit@copula
            
            ## simulate n_sims sequences from copula
            max_prediction_horizon <-
                last_analysis_time_season_week + 1 -
                analysis_time_season_week
            sim_sequences <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
            na_rows <- seq_len(n_sims)
#    while(length(na_rows) > 0) {
#            for(sim_ind in na_rows) {
#            ## generate random parameters from estimation distribution
#            ## Note that the variance estimate for parameters is low; at least we're
#            ## accounting for some uncertainty though...
#            predictive_copula@parameters <- rmvnorm(1, copula_fit@estimate, sigma = copula_fit@var.est)[1, ]
#            if(identical(class(predictive_copula)[1], "normalCopula")) {
#                ## randomly generated parameters above may not yield a positive definite correlation matrix; correct
#                predictive_copula@parameters <- makePosDef(getSigma(predictive_copula))[1, 1 + seq_along(predictive_copula@parameters)]
#            }
#                
#                ## simulate sequence from copula
#                sim_sequences[sim_ind, ] <- rCopula(1, predictive_copula)[1, ]
#            }
    
            
            sim_sequences <- rCopula(n_sims, predictive_copula)
            ## NA rows result when random parameters are problematic
#        na_rows <- which(apply(sim_sequences, 1, function(ss_row) any(is.na(ss_row))))
#    }
            
            ### get quantiles from marginal predictive distributions corresponding to
            ### values simulated from copula
            analysis_time_ind <- which(data$season == analysis_time_season &
                    data$season_week == analysis_time_season_week)
            trajectory_samples <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
            for(prediction_horizon in seq_len(max_prediction_horizon)) {
                trajectory_samples[, prediction_horizon] <-
                    kcde_predict(
                        p = sim_sequences[, prediction_horizon],
                        n = 100000,
                        kcde_fit = kcde_fits_by_prediction_horizon[[prediction_horizon]],
                        prediction_data =
                            data[seq_len(analysis_time_ind), , drop = FALSE],
                        leading_rows_to_drop = 0L,
                        trailing_rows_to_drop = 0L,
                        additional_training_rows_to_drop = NULL,
                        prediction_type = "quantile",
                        log = TRUE
                    )
                if(differencing) {
                    trajectory_samples[, prediction_horizon] <-
                        trajectory_samples[, prediction_horizon] *
                        data[analysis_time_ind + prediction_horizon - 52, prediction_target_var]
                }
            }
            
            ## Augment trajectory samples with previous observed incidence values
            season_start_ind <- which(data$season == analysis_time_season &
                    data$season_week == 1)
            if(season_start_ind < analysis_time_ind) {
                trajectory_samples <- cbind(
                    matrix(
                        rep(data[seq(from = season_start_ind, to = analysis_time_ind - 1), orig_prediction_target_var], each = n_sims),
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
            peak_week_height_by_sim_ind <- sapply(peak_week_height_by_sim_ind,
                function(height) {
                    which(incidence_bins$lower <= height &
                        incidence_bins$upper > height)
            })
            
            results[results_save_row, paste0("peak_height_", seq_len(n_sims))] <-
                peak_week_height_by_sim_ind
            
            ## Get log scores
            results[results_save_row, "peak_week_log_score"] <- log(sum(peak_week_by_sim_ind == observed_peak_week)) - log(n_sims)
            results[results_save_row, "peak_height_log_score"] <- log(sum(peak_week_height_by_sim_ind == observed_peak_height)) - log(n_sims)
        }
    }
    
    case_descriptor <- paste0(
        data_set,
        "-max_lag_", max_lag,
        "-max_seasonal_lag_", max_seasonal_lag,
        "-filtering_", filtering,
        "-differencing_", differencing,
        "-seasonality_", seasonality,
        "-bw_parameterization_", bw_parameterization
    )
    saveRDS(results,
        file.path(prediction_save_path,
            paste0("peak-week-", case_descriptor, ".rds")))
}

#library(plyr)
#library(dplyr)
#library(tidyr)
#library(ggplot2)
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
#    geom_line(aes(x = as.Date(prediction_time), y = sim_incidence, group = sim_ind), colour = "blue", alpha = 0.01, data = ts_for_plot) +
#    geom_hline(yintercept = 4.25906) +
#    geom_point(aes(x = as.Date(data$time[analysis_time_ind + peak_week_by_sim_ind - (analysis_time_season_week)]), y = peak_week_height_by_sim_ind)) +
#    scale_y_log10() +
#    theme_bw()
