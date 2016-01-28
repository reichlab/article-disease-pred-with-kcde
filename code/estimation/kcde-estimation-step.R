library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(doMC)

### Get command line arguments
args <- commandArgs(trailingOnly=TRUE)

## ILI, Dengue, or simulated data?
data_set <- args[1]

## Prediction horizon -- integer number of steps ahead to predict
prediction_horizon <- as.integer(args[2])

## Maximum number of lagged observations to explore using for prediction
max_lag <- as.integer(args[3])

## Include filtered cases as predictive variables?
filtering <- as.logical(args[4])

## Include terms capturing seasonality?
seasonality <- as.logical(args[5])

## Parameterization of bandwidth -- "diagonal" or "full"
bw_parameterization <- args[6]

## Path to save results in
save_path <- args[7]

## Manually set values for testing purposes
data_set <- "ili_national"
prediction_horizon <- 1L
max_lag <- 1L
filtering <- FALSE
seasonality <- FALSE
bw_parameterization <- "diagonal"
save_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-results"


### Load data set and set variables describing how the fit is performed
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
    continuous_var_names <- c("weighted_ili", "filtered_weighted_ili")
    discrete_var_names <- NULL
} else if(identical(data_set, "dengue_sj")) {
    ## Load data for Dengue fever in San Juan
} else if(identical(substr(data_set, 1, 3), "sim")) {
    ## Generate data for simulation study
    sim_ind <- as.integer(substr(data_set, 5, nchar(data_set)))
}



### Assemble control parameters for KCDE estimation process

## List describing kernel components -- depends on the values of
## prediction_horizon, filtering, seasonality, and bw_parameterization
kernel_components <- list()

## If requested, periodic kernel component capturing seasonality
if(seasonality) {
    kernel_components <- c(kernel_components,
        list(list(
            vars_and_offsets = data.frame(var_name = "time_index",
                offset_value = 0L,
                offset_type = "lag",
                combined_name = "time_index_lag0",
                stringsAsFactors = FALSE),
            kernel_fn = periodic_kernel,
            theta_fixed = list(period = pi / 365.2425), # 365.2425 is the mean number of days in a year
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_periodic_kernel,
            initialize_kernel_params_args = NULL,
            vectorize_kernel_params_fn = vectorize_params_periodic_kernel,
            vectorize_kernel_params_args = NULL,
            update_theta_from_vectorized_theta_est_fn =
                update_theta_from_vectorized_theta_est_periodic_kernel,
            update_theta_from_vectorized_theta_est_args = NULL
    )))
}

## Kernel components for observed values of incidence
## First step is setup: create list of data frames specifying groups of
## variables and offsets included in each kernel component
if(identical(bw_parameterization, "diagonal")) {
    ## Separate kernel components for each prediction target variable and
    ## predictive variable
    
    vars_and_offsets_groups <- list()
    
    ## Group of variable names and offsets for prediction target
    new_vars_and_offsets_group <- data.frame(
        var_name = prediction_target_var,
        offset_value = prediction_horizon,
        offset_type = "horizon",
        stringsAsFactors = FALSE
    )
    new_vars_and_offsets_group$combined_name <- paste0(
        new_vars_and_offsets_group$var_name,
        "_",
        new_vars_and_offsets_group$offset_type,
        new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- c(vars_and_offsets_groups,
        list(new_vars_and_offsets_group))
    
    ## Groups of variable names and offsets for lagged prediction target
    for(lag_value in seq(from = 0, to = max_lag)) {
        ## Group for lagged "raw"/unfiltered observed incidence
        new_vars_and_offsets_group <- data.frame(
            var_name = prediction_target_var,
            offset_value = lag_value,
            offset_type = "lag",
            stringsAsFactors = FALSE
        )
        new_vars_and_offsets_group$combined_name <- paste0(
            new_vars_and_offsets_group$var_name,
            "_",
            new_vars_and_offsets_group$offset_type,
            new_vars_and_offsets_group$offset_value
        )
        vars_and_offsets_groups <- c(vars_and_offsets_groups,
            list(new_vars_and_offsets_group))
        
        ## If requested, group for lagged filtered observed incidence
        if(filtering) {
            new_vars_and_offsets_group <- data.frame(
                var_name = paste0("filtered_", prediction_target_var),
                offset_value = lag_value,
                offset_type = "lag",
                stringsAsFactors = FALSE
            )
            new_vars_and_offsets_group$combined_name <- paste0(
                new_vars_and_offsets_group$var_name,
                "_",
                new_vars_and_offsets_group$offset_type,
                new_vars_and_offsets_group$offset_value
            )
            vars_and_offsets_groups <- c(vars_and_offsets_groups,
                list(new_vars_and_offsets_group))
        }
    }
} else if(identical(bw_parameterization, "full")) {
    ## One kernel component for prediction target variable and all predictive
    ## variables
    
    ## Prediction target variable
    new_vars_and_offsets_group <- data.frame(
        var_name = prediction_target_var,
        offset_value = prediction_horizon,
        offset_type = "horizon",
        stringsAsFactors = FALSE
    )
    
    ## Lagged prediction target == predictive variables
    for(lag_value in seq(from = 0, to = max_lag)) {
        ## Lagged "raw"/unfiltered observed incidence
        new_vars_and_offsets_group <- rbind(
            new_vars_and_offsets_group,
            data.frame(
                var_name = prediction_target_var,
                offset_value = lag_value,
                offset_type = "lag",
                stringsAsFactors = FALSE
            )
        )
        
        ## If requested, lagged filtered incidence
        if(filtering) {
            new_vars_and_offsets_group <- rbind(
                new_vars_and_offsets_group,
                data.frame(
                    var_name = paste0("filtered_", prediction_target_var),
                    offset_value = lag_value,
                    offset_type = "lag",
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    
    ## Add combined_name column and put in a list for further processing below
    new_vars_and_offsets_group$combined_name <- paste0(
        new_vars_and_offsets_group$var_name,
        "_",
        new_vars_and_offsets_group$offset_type,
        new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- list(new_vars_and_offsets_group)
} else {
    stop("Invalid bandwidth parameterization")
}

## Second step is to actually append the kernel component descriptions to the
## kernel_components list
## THIS CODE NEEDS TO BE UPDATED TO GET discrete_var_range_fns WHEN NEEDED
kernel_components <- c(kernel_components,
    lapply(vars_and_offsets_groups, function(vars_and_offsets) {
        lower_trunc_bds <- rep(-Inf, nrow(vars_and_offsets))
        names(lower_trunc_bds) <- vars_and_offsets$combined_name
        upper_trunc_bds <- rep(Inf, nrow(vars_and_offsets))
        names(upper_trunc_bds) <- vars_and_offsets$combined_name
        
        return(list(
            vars_and_offsets = vars_and_offsets,
            kernel_fn = log_pdtmvn_kernel,
            rkernel_fn = rlog_pdtmvn_kernel,
            theta_fixed = list(
                parameterization = "bw-chol-decomp",
                continuous_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$var_name %in% continuous_var_names],
                discrete_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$var_name %in% discrete_var_names],
                discrete_var_range_fns = NULL,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds
            ),
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
            initialize_kernel_params_args = NULL,
            vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
            vectorize_kernel_params_args = NULL,
            update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
            update_theta_from_vectorized_theta_est_args = NULL
        ))
    })
)


## Set up filter_control
if(filtering) {
    ## filter the prediction target variable
    filter_control <- list(
        list(
            var_name = prediction_target_var,
            max_filter_window_size = 13,
            filter_fn = two_pass_signal_filter,
            fixed_filter_params = list(
                n = 12,
                type = "pass",
                impute_fn = interior_linear_interpolation
            ),
            filter_args_fn = compute_filter_args_butterworth_filter,
            filter_args_args = NULL,
            initialize_filter_params_fn =
                initialize_filter_params_butterworth_filter,
            initialize_filter_params_args = NULL,
            vectorize_filter_params_fn =
                vectorize_filter_params_butterworth_filter,
            vectorize_filter_params_args = NULL,
            update_filter_params_from_vectorized_fn =
                update_filter_params_from_vectorized_butterworth_filter,
            update_filter_params_from_vectorized_args = NULL,
            transform_fn = log,
            detransform_fn = exp
        )
    )
} else {
    ## no filtering
    filter_control <- NULL
}


## Assemble kernel_components and filter_control created above,
## along with other parameters controling KCDE definition and estimation
kcde_control <- create_kcde_control(X_names = "time_index",
    y_names = prediction_target_var,
    time_name = "time",
    prediction_horizons = prediction_horizon,
    filter_control = filter_control,
    kernel_components = kernel_components,
    crossval_buffer = ymd("2010-01-01") - ymd("2009-01-01"),
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_args = list(
        prediction_type = "distribution",
        log = TRUE),
    loss_args = NULL)



### Do estimation
registerDoMC(cores=3)
fit_time <- system.time({
    kcde_fit <- kcde(data = data,
        kcde_control = kcde_control)
})



### Save results
saveRDS(kcde_fit,
    file = file.path(save_path,
        paste0("kcde_fit_",
            data_set,
            "_max_lag", as.character(max_lag),
            "_filtering", as.character(filtering),
            "_seasonality", as.character(seasonality),
            "_bw_parameterization", bw_parameterization,
            ".rds")
    )
)

saveRDS(fit_time,
    file = file.path(save_path,
        paste0("fit_time_",
            data_set,
            "_max_lag", as.character(max_lag),
            "_filtering", as.character(filtering),
            "_seasonality", as.character(seasonality),
            "_bw_parameterization", bw_parameterization,
            ".rds")
    )
)
