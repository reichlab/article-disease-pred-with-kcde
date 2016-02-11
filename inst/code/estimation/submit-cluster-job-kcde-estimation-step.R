options(warn = 2, error = recover)

all_data_sets <- c("ili_national")
all_prediction_horizons <- as.character(seq_len(52))
all_max_lags <- as.character(c(1L))
all_filtering_values <- c("FALSE", "TRUE")
all_seasonality_values <- c("FALSE", "TRUE")
all_bw_parameterizations <- c("diagonal", "full")

#all_prediction_horizons <- as.character(1L)
#all_max_lags <- as.character(c(1L))
#all_filtering_values <- c("TRUE")
#all_seasonality_values <- c("TRUE")
#all_bw_parameterizations <- c("full")

mem_req <- "1000"
time_req <- "60:00"
queue_req <- "long"
lsfoutfilename <- "kcde-est.out"

for(data_set in all_data_sets) {
    if(identical(data_set, "ili_national")) {
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
        results_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-results"
        scripts_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-scripts"
    }
    
    for(prediction_horizon in all_prediction_horizons) {
        for(max_lag in all_max_lags) {
            for(filtering in all_filtering_values) {
                for(seasonality in all_seasonality_values) {
                    cores_req <- as.character((as.numeric(max_lag) + 1) * (as.numeric(as.logical(filtering)) + 1) + as.numeric(as.logical(seasonality)))
                    
                    requestCmds <- "#!/bin/bash\n"
                    requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
                            "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
                            "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
                            "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
                            "#BSUB -W ", time_req, " # run time\n",
                            "#BSUB -q ", queue_req, " # which queue we want to run in\n")

                    for(bw_parameterization in all_bw_parameterizations) {
                        case_descriptor <- paste0(
                            data_set,
                            "-prediction_horizon_", prediction_horizon,
                            "-max_lag_", max_lag,
                            "-filtering_", filtering,
                            "-seasonality_", seasonality,
                            "-bw_parameterization_", bw_parameterization
                        )
                        
                        filename <- paste0(scripts_path, "/submit-kcde-estimation-step-", case_descriptor, ".sh")
                        
                        cat(requestCmds, file = filename)
                        cat("module load R/3.2.2\n", file = filename, append = TRUE)
                        cat(paste0("R CMD BATCH --vanilla \'--args ",
                            data_set, " ",
                            prediction_horizon, " ",
                            max_lag, " ",
                            filtering, " ",
                            seasonality, " ",
                            bw_parameterization, " ",
                            results_path,
                            "\'  /home/er71a/kcde-applied-paper/R/kcde-estimation-step.R ",
                            scripts_path, "/output-kde-estimation-step-", case_descriptor, ".Rout"),
                            file = filename, append = TRUE)
                        
                        bsubCmd <- paste0("bsub < ", filename)
                        
                        system(bsubCmd)
                    } # bw_parameterization
                } # seasonality
            } # filtering
        } # max_lag
    } # prediction_horizon
} # data_set
