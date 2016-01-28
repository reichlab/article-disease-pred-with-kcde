options(warn = 2, error = recover)

all_data_sets <- c("ili_national")
all_prediction_horizons <- as.character(seq_len(52))
all_max_lags <- as.character(c(3L))
all_filtering_values <- c("FALSE", "TRUE")
all_seasonality_values <- c("FALSE", "TRUE")
all_bw_parameterizations <- c("diagonal", "full")

all_prediction_horizons <- as.character(1L)
all_max_lags <- as.character(c(1L))
all_filtering_values <- c("FALSE")
all_seasonality_values <- c("FALSE")
all_bw_parameterizations <- c("diagonal")


for(data_set in all_data_sets) {
    if(identical(data_set, "ili_national")) {
        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
    }
    
    for(prediction_horizon in all_prediction_horizons) {
        for(max_lag in all_max_lags) {
            for(filtering in all_filtering_values) {
                for(seasonality in all_seasonality_values) {
                    for(bw_parameterization in all_bw_parameterizations) {
                        system_cmd <- paste0("R --vanilla --args ",
                            data_set, " ",
                            prediction_horizon, " ",
                            max_lag, " ",
                            filtering, " ",
                            seasonality, " ",
                            bw_parameterization, " ",
                            save_path,
                            " < /media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/kcde-estimation-step.R")
                        
                        system(system_cmd, intern = TRUE)
                    } # bw_parameterization
                } # seasonality
            } # filtering
        } # max_lag
    } # prediction_horizon
} # data_set
