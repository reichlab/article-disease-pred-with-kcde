library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(copula)
library(mvtnorm)

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


## set up
data_set <- "ili_national"
max_lag <- 1L
max_seasonal_lag <- 0L
filtering <- FALSE
differencing <- FALSE
seasonality <- TRUE
bw_parameterization <- "full"

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
        prediction_target_var <- "weighted_ili_ratio"
        train_seasons <- paste0(seq(from = 1998, to = 2009), "/", seq(from = 1999, to = 2010))
    } else {
        prediction_target_var <- "weighted_ili"
        train_seasons <- paste0(seq(from = 1997, to = 2009), "/", seq(from = 1998, to = 2010))
    }
    
    kernel_fn <- log_pdtmvn_kernel
    rkernel_fn <- rlog_pdtmvn_kernel
    
    variable_selection_method <- "all_included"
    crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")
    
    season_length <- 33L
    analysis_seasons <- c("2012/2013", "2013/2014", "2014/2015")
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
} else if(identical(data_set, "dengue_sj")) {
    ## Load data for Dengue fever in San Juan
    data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data/San_Juan_Training_Data.csv")
    
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
    
    kernel_fn <- log_pdtmvn_kernel
    rkernel_fn <- rlog_pdtmvn_kernel
    
    variable_selection_method <- "all_included"
    crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")
}

## generate peak week timing and height estimates
for(analysis_time_season in analysis_seasons) {
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
        for(sim_ind in na_rows) {
            ## generate random parameters from estimation distribution
            ## Note that the variance estimate for parameters is low; at least we're
            ## accounting for some uncertainty though...
            predictive_copula@parameters <- rmvnorm(1, copula_fit@estimate, sigma = copula_fit@var.est)[1, ]
            if(identical(class(predictive_copula)[1], "normalCopula")) {
                ## randomly generated parameters above may not yield a positive definite correlation matrix; correct
                predictive_copula@parameters <- makePosDef(getSigma(predictive_copula))[1, 1 + seq_along(predictive_copula@parameters)]
            }
            
            ## simulate sequence from copula
            sim_sequences[sim_ind, ] <- rCopula(1, predictive_copula)[1, ]
        }
        
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
}
}

ribbons_for_plot <- 

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

ts_for_plot <- trajectory_samples %>%
    as.data.frame() %>%
    `colnames<-`(as.character(seq_len(ncol(trajectory_samples)))) %>%
    mutate(sim_ind = seq_len(nrow(trajectory_samples))) %>%
    gather_("prediction_horizon", "sim_incidence", as.character(seq_len(ncol(trajectory_samples)))) %>%
    mutate(prediction_time = data$time[analysis_time_ind + as.integer(prediction_horizon)])

ggplot() +
    geom_line(aes(x = as.Date(time), y = weighted_ili), data = data[seq(from = analysis_time_ind - 10, to = analysis_time_ind + 52), ]) +
    geom_line(aes(x = as.Date(prediction_time), y = sim_incidence, group = sim_ind), alpha = 0.01, data = ts_for_plot) +
    theme_bw()

