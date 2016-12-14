library(surveillance)
library(lubridate)

## Note that the surveillance package applies only to discrete count data
## so we apply it only to the dengue data, not influenza
data_set <- "dengue_sj"

### Load data set and set variables describing how the fit is performed
## Load data for Dengue fever in San Juan
data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Testing_Data.csv")

## Restrict to data from 1990/1991 through 2008/2009 seasons
train_seasons <- paste0(1990:2008, "/", 1991:2009)
train_inds <- which(data$season %in% train_seasons)

## convert dates
data$time <- ymd(data$week_start_date)

## Add time_index column.  This is used for calculating the periodic kernel.
## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
## The origin is arbitrary.
data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))


## data in format for surveillance package
train_data <- sts(data$total_cases,
    start = c(year(data[1, "time"]), week(data[1, "time"])),
    freq = 52L)

## data frame with model specifications to evaluate
family_values <- c("Poisson", "NegBin1")
S_ar_values <- 0:3
S_end_values <- 0:3
lag_ar_values <- 1:3

model_specifications <- as.data.frame(
  expand.grid(
    family = family_values,
    S_ar = S_ar_values,
    S_end = S_end_values,
    lag_ar = lag_ar_values,
    mean_log_score = NA_real_,
    stringsAsFactors = FALSE),
  stringsAsFactors = FALSE)

fits <- vector("list", nrow(model_specifications))

## fit: note, no model component for neighbors since we're looking at a univariate time series
for(specification_ind in seq_len(nrow(model_specifications))) {
  family <- model_specifications$family[specification_ind]
  S_ar <- model_specifications$S_ar[specification_ind]
  lag_ar <- model_specifications$lag_ar[specification_ind]
  S_end <- model_specifications$S_end[specification_ind]
  
  fits[[specification_ind]] <- hhh4(train_data,
              control = list(
                ar = list(f = addSeason2formula(f = ~ 1, S = S_ar, period = 52), lag = lag_ar),
                end = list(f = addSeason2formula(f = ~ 1, S = S_end, period = 52)),
                subset = seq(from = lag_ar + 1, to = train_inds[length(train_inds)]),
                family = family
              ))
  
  ## evaluate via one-step-ahead predictions as discussed in Held and Paul
  one_step_ahead_preds <- oneStepAhead(fits[[specification_ind]],
    tp = nrow(data) - 2*52)
  pred_scores <- scores(one_step_ahead_preds)
  
  model_specifications$mean_log_score[specification_ind] <- mean(pred_scores[, "logs"])
  model_specifications$mean_rps[specification_ind] <- mean(pred_scores[, "rps"])
  model_specifications$mean_dss[specification_ind] <- mean(pred_scores[, "dss"])
  model_specifications$mean_ses[specification_ind] <- mean(pred_scores[, "ses"])
}


surveillance_fits <- list(
  model_specifications = model_specifications,
  model_fits = fits
)

saveRDS(surveillance_fits,
    file = file.path(
        "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
        data_set,
        "estimation-results/surveillance-fits.rds"))
