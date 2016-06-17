options(warn = 2, error = recover)

all_data_sets <- c("ili_national", "dengue_sj", "sim")

cores_req <- "4"
mem_req <- "1000"
time_req <- "100:00"
queue_req <- "long"

for(data_set in all_data_sets) {
    if(identical(data_set, "ili_national")) {
        all_prediction_horizons <- as.character(seq_len(52)) ## Obtain fits for each prediction horizon from 1 to 52 weeks
        all_max_lags <- as.character(c(1L)) # use incidence at times t^* and t^* - 1 to predict incidence after t^*
        all_max_seasonal_lags <- as.character(0L) # not used
        all_filtering_values <- c("FALSE") # not used
        all_differencing_values <- "FALSE" # not used
        all_seasonality_values <- c("FALSE", "TRUE") # specifications without and with periodic kernel
        all_bw_parameterizations <- c("diagonal", "full") # specifications with diagonal and full bandwidth
        all_sim_n <- "NA" # not used for applications
        all_sim_families <- "NA" # not used for applications
        all_sim_run_inds <- 1L # not used for applications
        
        results_path <- "/home/er71a/kcde-applied-paper/R/application-influenza-scott-rule-start/estimation-results"
        scripts_path <- "/home/er71a/kcde-applied-paper/R/application-influenza-scott-rule-start/estimation-scripts"
    } else if(identical(data_set, "dengue_sj")) {
        all_prediction_horizons <- as.character(seq_len(52)) ## Obtain fits for each prediction horizon from 1 to 52 weeks
        all_max_lags <- as.character(c(1L)) # use incidence at times t^* and t^* - 1 to predict incidence after t^*
        all_max_seasonal_lags <- as.character(0L) # not used
        all_filtering_values <- c("FALSE") # not used
        all_differencing_values <- "FALSE" # not used
        all_seasonality_values <- c("FALSE", "TRUE") # specifications without and with periodic kernel
        all_bw_parameterizations <- c("diagonal", "full") # specifications with diagonal and full bandwidth
        all_sim_n <- "NA" # not used for applications
        all_sim_families <- "NA" # not used for applications
        all_sim_run_inds <- 1L # not used for applications
        
        results_path <- "/home/er71a/kcde-applied-paper/R/application-dengue-scott-rule-start/estimation-results"
        scripts_path <- "/home/er71a/kcde-applied-paper/R/application-dengue-scott-rule-start/estimation-scripts"
    } else if(identical(data_set, "sim")) {
        all_prediction_horizons <- "0" # not used for simulation
        all_max_lags <- "0" # not used for simulation
        all_max_seasonal_lags <- "0" # not used for simulation
        all_filtering_values <- "FALSE" # not used for simulation
        all_differencing_values <- "FALSE" # not used for simulation
        all_seasonality_values <- "FALSE" # not used for simulation
        
        all_bw_parameterizations <- c("diagonal", "full") # specifications with diagonal and full bandwidth
        all_sim_n <- c("100", "1000") # sample size for training set
        all_sim_families <- "multivariate-2d-discretized" # distribution data are sampled from
        all_sim_run_inds <- seq(from = 1, to = 500) # index for simulation trial
    } else {
        stop("Invalid data set")
    }
    
    for(sim_family in all_sim_families) {
        if(identical(data_set, "sim")) {
            results_path <- paste0("/home/er71a/kcde-applied-paper/R/sim-", sim_family, "-scott-rule-start/estimation-results")
            scripts_path <- paste0("/home/er71a/kcde-applied-paper/R/sim-", sim_family, "-scott-rule-start/estimation-scripts")
        }
        for(sim_n in all_sim_n) {
            for(sim_run_ind in all_sim_run_inds) {
                for(prediction_horizon in all_prediction_horizons) {
                    for(max_lag in all_max_lags) {
                        for(max_seasonal_lag in all_max_seasonal_lags) {
                            for(filtering in all_filtering_values) {
                                for(differencing in all_differencing_values) {
                                    for(seasonality in all_seasonality_values) {
                                        for(bw_parameterization in all_bw_parameterizations) {
                                            if(identical(data_set, "ili_national")) {
                                                data_set_and_sim_run_ind <- data_set
                                                lsfoutfilename <- "kcde-est-applications-scott-rule-start.out"
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
                                            } else if(identical(data_set, "dengue_sj")) {
                                                data_set_and_sim_run_ind <- data_set
                                                lsfoutfilename <- "kcde-est-application-dengue-scott-rule-start.out"
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
                                            } else {
                                                data_set_and_sim_run_ind <- paste0(data_set, "_", sim_run_ind)
                                                lsfoutfilename <- "kcde-est-simstudies-scott-rule-start.out"
                                                case_descriptor <- paste0(
                                                    data_set_and_sim_run_ind,
                                                    "-prediction_horizon_", prediction_horizon,
                                                    "-max_lag_", max_lag,
                                                    "-filtering_", filtering,
                                                    "-seasonality_", seasonality,
                                                    "-bw_parameterization_", bw_parameterization,
                                                    "-sim_n_", sim_n
                                                )
                                            }
                                            
                                            filename <- paste0(scripts_path, "/submit-kcde-estimation-step-", case_descriptor, ".sh")
					                        
                                            requestCmds <- "#!/bin/bash\n"
                                            requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
                                                "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
                                                "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
                                                "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
                                                "#BSUB -W ", time_req, " # run time\n",
                                                "#BSUB -q ", queue_req, " # which queue we want to run in\n")
                                            
                                            cat(requestCmds, file = filename)
                                            cat("module load R/3.2.2\n", file = filename, append = TRUE)
                                            cat(paste0("R CMD BATCH --vanilla \'--args ",
                                                    data_set_and_sim_run_ind, " ",
                                                    prediction_horizon, " ",
                                                    max_lag, " ",
                                                    max_seasonal_lag, " ",
                                                    filtering, " ",
                                                    differencing, " ",
                                                    seasonality, " ",
                                                    bw_parameterization, " ",
                                                    results_path, " ",
                                                    sim_n, " ",
                                                    sim_family,
                                                    "\'  /home/er71a/kcde-applied-paper/R/kcde-estimation-step.R ",
                                                    scripts_path, "/output-kde-estimation-step-", case_descriptor, ".Rout"),
                                                file = filename, append = TRUE)
                                            
                                            bsubCmd <- paste0("bsub < ", filename)
                                            
                                            system(bsubCmd)
                                        } # bw_parameterization
                                    } # seasonality
                                } # differencing
                            } # filtering
                        } # max_seasonal_lag
                    } # max_lag
                } # prediction_horizon
            } # sim_run_ind
        } # sim_n
    } # sim_family
} # data_set
