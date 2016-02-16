options(warn = 2, error = recover)

#all_data_sets <- c("ili_national")
all_data_sets <- c("sim")

mem_req <- "1000"
time_req <- "60:00"
queue_req <- "long"

for(data_set in all_data_sets) {
    if(identical(data_set, "ili_national")) {
        all_prediction_horizons <- as.character(seq_len(52))
        all_max_lags <- as.character(c(1L))
        all_filtering_values <- c("FALSE", "TRUE")
        all_seasonality_values <- c("FALSE", "TRUE")
        all_bw_parameterizations <- c("diagonal", "full")
        all_sim_n <- "NA"
        all_sim_families <- "NA"
        all_sim_run_inds <- 1L
        
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
        results_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-results"
        scripts_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-scripts"
    } else if(identical(data_set, "dengue_sj")) {
        all_prediction_horizons <- as.character(seq_len(52))
        all_max_lags <- as.character(c(1L))
        all_filtering_values <- c("FALSE", "TRUE")
        all_seasonality_values <- c("FALSE", "TRUE")
        all_bw_parameterizations <- c("diagonal", "full")
        all_sim_n <- "NA"
        all_sim_families <- "NA"
        all_sim_run_inds <- 1L
        
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
        results_path <- "/home/er71a/kcde-applied-paper/R/application-dengue/estimation-results"
        scripts_path <- "/home/er71a/kcde-applied-paper/R/application-dengue/estimation-scripts"
    } else if(identical(data_set, "sim")) {
        all_prediction_horizons <- "0"
        all_max_lags <- "0"
        all_filtering_values <- "FALSE"
        all_seasonality_values <- "FALSE"
        all_bw_parameterizations <- c("diagonal", "full")
        all_sim_n <- c("100", "1000")
        all_sim_families <- c("bivariate-B-discretized", "bivariate-C-discretized", "multivariate-2d-discretized", "multivariate-4d-discretized", "multivariate-6d-discretized")
        all_sim_run_inds <- seq(from = 2, to = 100)
    } else {
        stop("Invalid data set")
    }
    
    for(sim_family in all_sim_families) {
        if(identical(data_set, "sim")) {
#        save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/R/application-influenza/estimation-results"
            results_path <- paste0("/home/er71a/kcde-applied-paper/R/sim-", sim_family, "/estimation-results")
            scripts_path <- paste0("/home/er71a/kcde-applied-paper/R/sim-", sim_family, "/estimation-scripts")
        }
        for(sim_n in all_sim_n) {
            for(sim_run_ind in all_sim_run_inds) {
                for(prediction_horizon in all_prediction_horizons) {
                    for(max_lag in all_max_lags) {
                        for(filtering in all_filtering_values) {
                            for(seasonality in all_seasonality_values) {
                                for(bw_parameterization in all_bw_parameterizations) {
                                    if(data_set %in% c("ili_national", "dengue_sj")) {
                                        data_set_and_sim_run_ind <- data_set
                                        lsfoutfilename <- "kcde-est-applications.out"
                                        case_descriptor <- paste0(
                                            data_set,
                                            "-prediction_horizon_", prediction_horizon,
                                            "-max_lag_", max_lag,
                                            "-filtering_", filtering,
                                            "-seasonality_", seasonality,
                                            "-bw_parameterization_", bw_parameterization
                                        )
                                    } else {
                                        data_set_and_sim_run_ind <- paste0(data_set, "_", sim_run_ind)
                                        lsfoutfilename <- "kcde-est-simstudies.out"
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
                                    
                                    if(!file.exists(filename)) {
                                        cores_req <- as.character((as.numeric(max_lag) + 1) * (as.numeric(as.logical(filtering)) + 1) + as.numeric(as.logical(seasonality)))
                                        
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
                                                filtering, " ",
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
				    } # check if file exists
                                } # bw_parameterization
                            } # seasonality
                        } # filtering
                    } # max_lag
                } # prediction_horizon
            } # sim_run_ind
        } # sim_n
    } # sim_family
} # data_set
