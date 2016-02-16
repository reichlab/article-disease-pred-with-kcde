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
#data_set <- "ili_national"
#prediction_horizon <- 25L
#max_lag <- 1L
#filtering <- TRUE
#seasonality <- TRUE
#bw_parameterization <- "full"
#save_path <- "/home/er71a/kcde-applied-paper/R/application-influenza/estimation-results"
#data_set <- "dengue_sj"

data_set <- "sim_10"
prediction_horizon <- 0L
max_lag <- 0L
filtering <- FALSE
seasonality <- FALSE
bw_parameterization <- "full"
save_path <- "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/results/dengue_sj/estimation-results"
args <- c(rep("", 7), "100", "bivariate-B-discretized")
#args <- c(rep("", 7), "1000", "bivariate-C-discretized")
#args <- c(rep("", 7), "1000", "multivariate-4d-discretized")


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
    continuous_var_names <- c(
        paste0(c("weighted_ili", "filtered_weighted_ili"), "_horizon", rep(1:52, each=2)),
        paste0(c("weighted_ili", "filtered_weighted_ili"), "_lag", rep(seq(from = 0, to = max_lag), each=2))
    )
    discrete_var_names <- NULL
    predictive_vars <- paste0(c("weighted_ili", "filtered_weighted_ili"), "_lag", rep(seq(from = 0, to = max_lag), each=2))
    time_var <- "time"
    
    kernel_fn <- log_pdtmvn_kernel
    rkernel_fn <- rlog_pdtmvn_kernel
    
    variable_selection_method <- "stepwise"
    crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")
} else if(identical(data_set, "dengue_sj")) {
    ## Load data for Dengue fever in San Juan
    data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data/San_Juan_Training_Data.csv")
    
    ## Form variable with total cases + 1 which can be logged
    data$total_cases_plus_1 <- data$total_cases + 1
    
    ## convert dates
    data$time <- ymd(data$week_start_date)
    
    ## Add time_index column.  This is used for calculating the periodic kernel.
    ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
    ## The origin is arbitrary.
    data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
    
    prediction_target_var <- "total_cases_plus_1"
    continuous_var_names <- c(
        paste0(c("total_cases_plus_1", "filtered_total_cases_plus_1"), "_lag", rep(seq(from = 0, to = max_lag), each=2))
    )
    discrete_var_names <- c(
        paste0(c("total_cases_plus_1", "filtered_total_cases_plus_1"), "_horizon", rep(1:52, each=2))
    )
    predictive_vars <- paste0(c("total_cases_plus_1", "filtered_total_cases_plus_1"), "_lag", rep(seq(from = 0, to = max_lag), each=2))
    time_var <- "time"
    
    kernel_fn <- log_pdtmvn_kernel
    rkernel_fn <- rlog_pdtmvn_kernel
    
    variable_selection_method <- "stepwise"
    crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")
} else if(identical(substr(data_set, 1, 3), "sim")) {
    ## Load functions for generating simulated data
    source("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/code/sim-densities-sim-study-discretized-Duong-Hazelton.R")
#    source("/home/er71a/kcde-applied-paper/R/sim-densities-sim-study-discretized-Duong-Hazelton.R")
    
    ## Get arguments determining size of simulation and simulation family
    sim_n <- as.integer(args[8])
    sim_family <- args[9]
    sim_ind <- as.integer(substr(data_set, 5, nchar(data_set)))
    cat(sim_n)
    cat(sim_family)
    cat(sim_ind)
    cat("\n\n")
    cat(data_set)
    cat("\n")
    cat(substr(data_set, 5, nchar(data_set)))
    cat("\n")
    
    ## Set up random number generation
    ## We use the rstream package for RNG, with
    ## a separate stream of random numbers used for each combination of
    ## simulation index, simulation family, and observed sample size
    ## We allow for up to 10000 simulations for each combinatio nof
    ## simulation family and sample size
    rng_ind <- sim_ind
    cat(rng_ind)
    rng_ind <- rng_ind + switch(sim_family,
        "bivariate-A-discretized" = 0,
        "bivariate-B-discretized" = 20000,
        "bivariate-C-discretized" = 40000,
        "bivariate-D-discretized" = 60000,
        "multivariate-2d-discretized" = 80000,
        "multivariate-4d-discretized" = 100000,
        "multivariate-6d-discretized" = 120000,
        stop("Invalid simulation family.")
    )
    cat(rng_ind)
    rng_ind <- rng_ind + switch(as.character(sim_n),
        "100" = 0,
        "1000" = 10000,
        stop("Invalid simulation sample size.")
    )
    cat(rng_ind)
    
    set.seed(51158) # this was randomly generated
    
    library(rstream)
    rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))
    for(i in seq_len(rng_ind)) {
        rstream.nextsubstream(rngstream)
    }
    
    rstream.RNG(rngstream)
    
    ## Generate data
    data <- sim_from_pdtmvn_mixt(n = sim_n, sim_family = sim_family)
    
    prediction_target_var <- paste0("X", ncol(data))
    predictive_vars <- paste0("X", seq_len(ncol(data) - 1))
    
    continuous_var_names <- NULL
    discrete_var_names <- paste0(colnames(data), c(rep("_lag", ncol(data) - 1), "_horizon"), 0)
    
    time_var <- NULL
    
    kernel_fn <- pdtmvn_kernel
    rkernel_fn <- rpdtmvn_kernel
    
    variable_selection_method <- "all_included"
    crossval_buffer <- 0L
} else {
    stop("Invalid data set argument.")
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
    
    ## Groups of variable names and offsets for lagged predictive variables
    
    for(lag_value in seq(from = 0, to = max_lag)) {
        for(predictive_var in predictive_vars) {
            ## Group for lagged "raw"/unfiltered observed incidence
            new_vars_and_offsets_group <- data.frame(
                var_name = predictive_var,
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
                    var_name = paste0("filtered_", predictive_var),
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
        for(predictive_var in predictive_vars) {
            ## Lagged "raw"/unfiltered observed incidence
            new_vars_and_offsets_group <- rbind(
                new_vars_and_offsets_group,
                data.frame(
                    var_name = predictive_var,
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
                        var_name = paste0("filtered_", predictive_var),
                        offset_value = lag_value,
                        offset_type = "lag",
                        stringsAsFactors = FALSE
                    )
                )
            }
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

if(data_set %in% c("ili_national", "dengue_sj")) {
    #' Compute whether each element of x is equal to the log of an integer,
    #' up to a specified tolerance level.
    #' 
    #' @param x numeric
    #' @param tolerance numeric tolerance for comparison of integer values
    #' 
    #' @return logical vector of same length as x; entry i is TRUE if
    #'     x[i] is within tol of as.integer(x[i])
    equals_log_integer <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
        return(sapply(x, function(x_i) {
                    return(isTRUE(all.equal(x_i, log(as.integer(exp(x_i))))))
                }))
    }
    
    #' Compute log(exp(x) - 0.5)
    #' Used as default "a" function
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_exp_x_minus_0.5 <- function(x) {
        return(log(exp(x) - 0.5))
    }
    
    #' Compute log(exp(x) + 0.5)
    #' Used as default "a" function
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_exp_x_plus_0.5 <- function(x) {
        return(log(exp(x) + 0.5))
    }
    
    #' Compute log(round(exp(x))) in such a way that the rounding function always rounds up _.5
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_round_up_0.5_exp <- function(x) {
        exp_x <- exp(x)
        
        inds_ceil <- exp_x - floor(exp_x) >= 0.5
        
        exp_x[inds_ceil] <- ceiling(exp_x[inds_ceil])
        exp_x[!inds_ceil] <- floor(exp_x[!inds_ceil])
        
        return(log(exp_x))
    }
    
    
    a_fn <- log_exp_x_minus_0.5
    b_fn <- log_exp_x_plus_0.5
    in_range_fn <- equals_log_integer
    discretizer_fn <- log_round_up_0.5_exp
} else {
    a_fn <- x_minus_0.25
    b_fn <- x_plus_0.25
    in_range_fn <- equals_half_integer
    discretizer_fn <- round_to_half_integer
}

kernel_components <- c(kernel_components,
    lapply(vars_and_offsets_groups, function(vars_and_offsets) {
        lower_trunc_bds <- rep(-Inf, nrow(vars_and_offsets))
        names(lower_trunc_bds) <- vars_and_offsets$combined_name
        upper_trunc_bds <- rep(Inf, nrow(vars_and_offsets))
        names(upper_trunc_bds) <- vars_and_offsets$combined_name
        
        if("horizon" %in% vars_and_offsets$offset_type || identical(substr(data_set, 1, 3), "sim")) {
            discrete_var_range_fns <- lapply(
                 discrete_var_names,
                 function(discrete_var_name) {
                     list(a = a_fn,
                         b = b_fn,
                         in_range = in_range_fn,
                         discretizer = discretizer_fn)
                 }
            ) %>%
                `names<-`(discrete_var_names)
        } else {
            discrete_var_range_fns <- NULL
        }
        
        return(list(
            vars_and_offsets = vars_and_offsets,
            kernel_fn = kernel_fn,
            rkernel_fn = rkernel_fn,
            theta_fixed = list(
                parameterization = "bw-chol-decomp",
                continuous_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$combined_name %in% continuous_var_names],
                discrete_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$combined_name %in% discrete_var_names],
                discrete_var_range_fns = discrete_var_range_fns,
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
    time_name = time_var,
    prediction_horizons = prediction_horizon,
    filter_control = filter_control,
    kernel_components = kernel_components,
    crossval_buffer = crossval_buffer,
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_args = list(
        prediction_type = "distribution",
        log = TRUE),
    loss_args = NULL,
    variable_selection_method = variable_selection_method)



### Do estimation
num_cores <- (max_lag + 1) * (filtering + 1) + seasonality
registerDoMC(cores = num_cores)
fit_time <- system.time({
    kcde_fit <- kcde(data = data,
        kcde_control = kcde_control)
})



### Save results
case_descriptor <- paste0(
    data_set,
    "-prediction_horizon_", prediction_horizon,
    "-max_lag_", max_lag,
    "-filtering_", filtering,
    "-seasonality_", seasonality,
    "-bw_parameterization_", bw_parameterization
)

if(identical(substr(data_set, 1, 3), "sim")) {
    case_descriptor <- paste0(
        case_descriptor,
        "-sim_n_", sim_n,
        "-sim_ind_", sim_ind)
}

saveRDS(kcde_fit,
    file = file.path(save_path,
        paste0("kcde_fit-",
            case_descriptor,
            ".rds")
    )
)

saveRDS(fit_time,
    file = file.path(save_path,
        paste0("fit_time-",
            case_descriptor,
            ".rds")
    )
)
