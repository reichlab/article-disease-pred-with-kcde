library(rstan)
library(plyr)
library(dplyr)
library(tidyr)
library(splines)
library(ggplot2)

ili_prediction_results_sarima <- readRDS("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/ili_national/prediction-results/sarima-predictions.rds")
ili_prediction_results_kcde <- readRDS("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/ili_national/prediction-results/kcde-predictions.rds")
ili_prediction_results <- rbind.fill(ili_prediction_results_sarima[!is.na(ili_prediction_results_sarima$log_score), ],
    ili_prediction_results_kcde)
ili_prediction_results$AE <- unlist(ili_prediction_results$AE)

ili_prediction_results$full_model_descriptor <- paste0(ili_prediction_results$model,
    "-seasonal_lag_", ili_prediction_results$max_seasonal_lag,
    "-differencing_", ili_prediction_results$differencing,
    "-periodic_", ili_prediction_results$seasonality,
    "-bw_", ili_prediction_results$bw_parameterization)


ili_prediction_log_score_diffs_wide <- ili_prediction_results %>%
    select(full_model_descriptor, prediction_time, prediction_horizon, log_score) %>%
    spread(full_model_descriptor, log_score)

ili_prediction_log_score_diffs_wide[, unique(ili_prediction_results$full_model_descriptor)] <-
    ili_prediction_log_score_diffs_wide[, unique(ili_prediction_results$full_model_descriptor)] -
    ili_prediction_log_score_diffs_wide[, "SARIMA-seasonal_lag_NA-differencing_NA-periodic_NA-bw_NA"]

ili_prediction_log_score_diffs_long <- ili_prediction_log_score_diffs_wide %>%
    gather_("model", "log_score_difference", unique(ili_prediction_results$full_model_descriptor)) %>%
    filter(model != "SARIMA-seasonal_lag_NA-differencing_NA-periodic_NA-bw_NA")

ili_prediction_log_score_diffs_long$model_indicator <- as.integer(as.factor(ili_prediction_log_score_diffs_long$model))

prediction_horizon_spline_basis_for_plot <- bs(1:52, df = 26, intercept = TRUE)
prediction_horizon_spline_basis_plot_df <- prediction_horizon_spline_basis_for_plot %>%
    as.data.frame() %>%
    mutate(prediction_horizon = 1:52) %>%
    gather_("spline", "value", as.character(1:26))

ggplot(prediction_horizon_spline_basis_plot_df) +
    geom_line(aes(x = prediction_horizon, y = value, colour = spline)) +
    theme_bw()



## Quantities to be passed to stan as data
N <- nrow(ili_prediction_log_score_diffs_long)
M <- 16
B_spline <- 26
B <- B_spline + 1

## Create spline basis
prediction_horizon_spline_basis <- cbind(
    matrix(1, nrow = N),
    bs(ili_prediction_log_score_diffs_long$prediction_horizon,
        df = B_spline,
        intercept = TRUE)
)

colnames(prediction_horizon_spline_basis) <- paste0("spline_basis_", seq_len(B))

## Spread spline basis across rows of an M * B matrix with basis functions for each model in a separate group of B columns
## There's got to be a slick way of doing this without a for loop, but I can't think of it
X <- matrix(0, nrow = N, ncol = M * B)
for(i in seq_len(N)) {
    X[i, seq_len(B) + (ili_prediction_log_score_diffs_long$model_indicator[i] - 1) * B] <- prediction_horizon_spline_basis[i, ]
}



prediction_time <- as.numeric(ili_prediction_log_score_diffs_long$prediction_time)
prediction_horizon <- ili_prediction_log_score_diffs_long$prediction_horizon
model <- ili_prediction_log_score_diffs_long$model_indicator
log_score_difference <- ili_prediction_log_score_diffs_long$log_score_difference



#X <- kronecker(diag(M), prediction_horizon_spline_basis)

rstan_options(auto_write = TRUE)
options(mc.cores = 2)

ili_results_fit <- stan(
    file = "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/code/postprocessing/influenza-results-model/influenza-results-model-spline-indep-stan.stan",
    data = c(
        "N",
        "M",
        "B",
        "X",
        "prediction_time",
#        "prediction_horizon",
#        "model",
        "log_score_difference"
    ),
    pars = c("beta", "sigma", "y_pred"),
    iter = 1000,
    chains = 2
)


saveRDS(ili_results_fit,
    file = "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/ili_national/analysis-results/influenza-results-model-spline-indep-fit-chains_2-iter_1000.rds")

ili_results_fit_sso <- as.shinystan(ili_results_fit)
launch_shinystan(ili_results_fit_sso)

ili_results_fit_samples <- as.data.frame(ili_results_fit)
saveRDS(ili_results_fit_samples,
    file = "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results/ili_national/analysis-results/influenza-results-model-samples-spline-indep-fit-chains_2-iter_1000.rds")


