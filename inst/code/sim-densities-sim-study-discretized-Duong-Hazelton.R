#' Compute whether each element of x is equal to an integer, or an integer +/- 0.5
#' up to a specified tolerance level.
#' 
#' @param x numeric
#' @param tolerance numeric tolerance for comparison of integer values
#' 
#' @return logical vector of same length as x; entry i is TRUE if
#'     x[i] is within tolerance of as.integer(x[i])
equals_half_integer <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
    return(sapply(x, function(x_i) {
        return(isTRUE(all.equal(x_i, as.integer(x_i), tolerance = tolerance)) ||
            isTRUE(all.equal(x_i, as.integer(x_i) - 0.5, tolerance = tolerance)) ||
            isTRUE(all.equal(x_i, as.integer(x_i) + 0.5, tolerance = tolerance))
        )
    }))
}

#' Compute x - 0.25
#' Used as "a" function
#' 
#' @param x numeric
#' 
#' @return x - 0.25
x_minus_0.25 <- function(x) {
    return(x - 0.25)
}

#' Compute x + 0.25
#' Used as "b" function
#' 
#' @param x numeric
#' 
#' @return x + 0.25
x_plus_0.25 <- function(x) {
    return(x + 0.25)
}

#' Round x to the nearest half integer in such a way that we always round up _.25
#' and _.75
#' 
#' @param x numeric
#' 
#' @return x rounded to the nearest half integer
round_to_half_integer <- function(x) {
    ## Get candidate return values
    x_int <- as.integer(x)
    round_candidates <- t(apply(as.matrix(x_int), 1, function(x_int_i) {
        return(x_int_i + seq(from = 1, to = -1, by = -0.5))
    }))
    
    ## Select column index in round_candidates with smallest difference from
    ## corresponding value of x
    abs_diff_from_round_candidates <- abs(sweep(round_candidates, 1, x))
    inds_keep <- apply(abs_diff_from_round_candidates, 1, which.min)
    
    return(round_candidates[cbind(seq_along(x), inds_keep)])
}

#' Get parameters for distributions used in simulation study
#' based on discretizing and conditioning on some covariates in the
#' simulation study from Duong and Hazelton (2005).  Parameters returned
#' are suitable for a call to kcde::simulate_values_from_product_kernel
#' 
#' @param sim_family a string describing the distribution to sample from:
#'   either "bivariate-A" through "bivariate-D" or
#'   "multivariate-2", "multivariate-4", or "multivariate-6"
#' 
#' @return list of parameters
get_dist_component_params_for_sim_family <- function(sim_family) {
    dist_component_params <- list()
    
    if(identical(sim_family, "bivariate-A-discretized")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(0.25, 0, 0, 1), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(2), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- 1
        
        dist_component_params$centers <- matrix(rep(0, 2), nrow = 1, ncol = 2) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
            kernel_fn = pdtmvn_kernel,
            rkernel_fn = rpdtmvn_kernel,
            theta_fixed = list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds
            ),
            vars_and_offsets = data.frame(combined_name = x_names)
        ))
        
        dist_component_params$theta <- list(list(
            parameterization = "bw-chol-decomp",
            continuous_vars = NULL,
            discrete_vars = x_names,
            discrete_var_range_fns = discrete_var_range_fns,
            lower = lower_trunc_bds,
            upper = upper_trunc_bds,
            x_names = x_names,
            bw = sigma,
            continuous_var_col_inds = NULL,
            discrete_var_col_inds = seq_along(x_names)
        ))
    } else if(identical(sim_family, "bivariate-A")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(0.25, 0, 0, 1), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- 1
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$centers <- matrix(rep(0, 2), nrow = 1, ncol = 2) %>%
            `colnames<-`(x_names)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else if(identical(sim_family, "bivariate-B-discretized")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(4/9, 0, 0, 4/9), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(2), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
            matrix(c(1, 0), nrow = 1, ncol = 2, byrow = TRUE),
            matrix(c(-1, 0), nrow = 1, ncol = 2, byrow = TRUE)
        ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "bivariate-B")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(4/9, 0, 0, 4/9), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
                matrix(c(1, 0), nrow = 1, ncol = 2, byrow = TRUE),
                matrix(c(-1, 0), nrow = 1, ncol = 2, byrow = TRUE)
            ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else if(identical(sim_family, "bivariate-C-discretized")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(1, 0.9, 0.9, 1), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(2), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
                matrix(c(1, -0.9), nrow = 1, ncol = 2, byrow = TRUE),
                matrix(c(-1, 0.9), nrow = 1, ncol = 2, byrow = TRUE)
            ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "bivariate-C")) {
        x_names <- paste0("X", seq_len(2))
        sigma <- matrix(c(1, 0.9, 0.9, 1), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
                matrix(c(1, -0.9), nrow = 1, ncol = 2, byrow = TRUE),
                matrix(c(-1, 0.9), nrow = 1, ncol = 2, byrow = TRUE)
            ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else if(identical(sim_family, "bivariate-D-discretized")) {
        stop("Invalid sim_family") # The parameters given in the text of Duong and Hazelton are incorrect.
        x_names <- paste0("X", seq_len(2))
        sigma1 <- matrix(c(25/64, 1/5, 1/5, 25/64), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        sigma2 <- matrix(c(25/64, -1/4, -1/4, 25/64), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        sigma3 <- matrix(c(15/32, -1 / (4 * sqrt(3)), -1 / (4 * sqrt(3)), 5/8), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(2), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
                matrix(c(73/64, -5/6), nrow = 1, ncol = 2, byrow = TRUE),
                matrix(c(7/32, -5/3), nrow = 1, ncol = 2, byrow = TRUE),
                matrix(c(87/64, -5/6), nrow = 1, ncol = 2, byrow = TRUE)
            ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma1, sigma2, sigma3)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma1,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "bivariate-D")) {
        stop("Invalid sim_family") # The parameters given in the text of Duong and Hazelton are incorrect.
        x_names <- paste0("X", seq_len(2))
        sigma1 <- matrix(c(25/64, 1/5, 1/5, 25/64), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        sigma2 <- matrix(c(25/64, -1/4, -1/4, 25/64), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        sigma3 <- matrix(c(15/32, -1 / (4 * sqrt(3)), -1 / (4 * sqrt(3)), 5/8), nrow = 2, ncol = 2) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- rep(0.5, 2)
        
        dist_component_params$centers <- rbind(
            matrix(c(73/64, -5/6), nrow = 1, ncol = 2, byrow = TRUE),
            matrix(c(7/32, -5/3), nrow = 1, ncol = 2, byrow = TRUE),
            matrix(c(87/64, -5/6), nrow = 1, ncol = 2, byrow = TRUE)
        ) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma1, sigma2, sigma3)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma1,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
        
    } else if(identical(sim_family, "multivariate-2d-discretized")) {
        dimension <- 2L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(dimension), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- 1
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "multivariate-2d")) {
        dimension <- 2L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- 1
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else if(identical(sim_family, "multivariate-4d-discretized")) {
        dimension <- 4L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(dimension), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- 1
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "multivariate-4d")) {
        dimension <- 4L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- 1
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else if(identical(sim_family, "multivariate-6d-discretized")) {
        dimension <- 6L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- lapply(seq_len(dimension), function(var_ind) {
                    list(a = x_minus_0.25,
                        b = x_plus_0.25,
                        in_range = equals_half_integer,
                        discretizer = round_to_half_integer
                    )
                }) %>%
            `names<-`(x_names)
        
        dist_component_params$weights <- 1
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = NULL,
                    discrete_vars = x_names,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = NULL,
                discrete_vars = x_names,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = NULL,
                discrete_var_col_inds = seq_along(x_names)
            ))
    } else if(identical(sim_family, "multivariate-6d")) {
        dimension <- 6L
        x_names <- paste0("X", seq_len(dimension))
        sigma <- matrix(0.9, nrow = dimension, ncol = dimension) %>%
            `rownames<-`(x_names) %>%
            `colnames<-`(x_names)
        diag(sigma) <- 1
        lower_trunc_bds <- rep(-Inf, length(x_names)) %>%
            `names<-`(x_names)
        upper_trunc_bds <- rep(Inf, length(x_names)) %>%
            `names<-`(x_names)
        
        discrete_var_range_fns <- NULL
        
        dist_component_params$weights <- 1
        
        dist_component_params$bws <- list(sigma, sigma)
        
        dist_component_params$centers <- matrix(rep(0, dimension), nrow = 1, ncol = dimension) %>%
            `colnames<-`(x_names)
        
        dist_component_params$kernel_components <- list(list(
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = list(
                    parameterization = "bw-chol-decomp",
                    continuous_vars = x_names,
                    discrete_vars = NULL,
                    discrete_var_range_fns = discrete_var_range_fns,
                    lower = lower_trunc_bds,
                    upper = upper_trunc_bds
                ),
                vars_and_offsets = data.frame(combined_name = x_names)
            ))
        
        dist_component_params$theta <- list(list(
                parameterization = "bw-chol-decomp",
                continuous_vars = x_names,
                discrete_vars = NULL,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                x_names = x_names,
                bw = sigma,
                continuous_var_col_inds = seq_along(x_names),
                discrete_var_col_inds = NULL
            ))
    } else {
        stop("Invalid sim_family")
    }
    
    return(dist_component_params)
}



sim_from_pdtmvn_mixt <- function(n, sim_family) {
    dist_component_params <- get_dist_component_params_for_sim_family(sim_family)
    
    result <- matrix(NA, nrow = n, ncol = ncol(dist_component_params$theta[[1]]$bw))
    colnames(result) <- paste0("X", seq_len(ncol(result)))
    
    sampled_kernel_inds <- sample(length(dist_component_params$weights),
        size = n,
        replace = TRUE,
        prob = dist_component_params$weights)
    
    for(kernel_ind in unique(sampled_kernel_inds)) {
        result_inds <- which(sampled_kernel_inds == kernel_ind)
        
        theta_for_sim <- dist_component_params$theta
        theta_for_sim$bw <- dist_component_params$bws[[kernel_ind]]
        result[result_inds, ] <- simulate_values_from_product_kernel(n = length(result_inds),
            center = dist_component_params$centers[kernel_ind, , drop = FALSE],
            kernel_components = dist_component_params$kernel_components,
            theta = theta_for_sim)[, colnames(dist_component_params$centers)]
    }
    
    return(result)
}



d_pdtmvn_mixt_conditional <- function(X, sim_family, conditional = TRUE, log = FALSE) {
    dist_component_params <- get_dist_component_params_for_sim_family(sim_family)
    
    log_kernel_component_values <- matrix(NA, nrow = nrow(X), ncol = length(dist_component_params$weights))
    
    kernel_fn_args <- dist_component_params$theta[[1]]
    for(ind in seq_len(ncol(log_kernel_component_values))) {
        kernel_fn_args$bw <- dist_component_params$bws[[ind]]
        kernel_fn_args$x <- dist_component_params$centers[ind, , drop = FALSE]
        kernel_fn_args$center <- X
        kernel_fn_args$log <- TRUE
        
        log_kernel_component_values[, ind] <-
            do.call(dist_component_params$kernel_components[[1]]$kernel_fn,
                kernel_fn_args)
        
        if(conditional) {
            conditioning_var_inds <- seq_len(ncol(kernel_fn_args$bw) - 1)
        
#            kernel_fn_args$bw <- dist_component_params$bws[[ind]][conditioning_var_inds, conditioning_var_inds, drop = FALSE]
            kernel_fn_args$x <- dist_component_params$centers[ind, conditioning_var_inds, drop = FALSE]
            kernel_fn_args$center <- X[, conditioning_var_inds, drop = FALSE]
            
            log_kernel_component_values[, ind] <- log_kernel_component_values[, ind] -
                do.call(dist_component_params$kernel_components[[1]]$kernel_fn,
                    kernel_fn_args)
        }
    }
    
    if(log) {
        return(apply(log_kernel_component_values, 1, sum))
    } else {
        return(exp(apply(log_kernel_component_values, 1, sum)))
    }
}
