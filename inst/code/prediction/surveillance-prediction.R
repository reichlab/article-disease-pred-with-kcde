library(surveillance)
library(cdcfluview)
library(plyr)
library(dplyr)
library(kcde)
library(lubridate)

selection_criteria <- "aic" # how to pick "best" hhh4 model
data_set <- "dengue_sj"
all_prediction_horizons <- 1:52
all_prediction_statistics <- c("log_score",
    "pt_pred",
    "AE",
    "interval_pred_lb_95",
    "interval_pred_ub_95",
    "interval_pred_lb_50",
    "interval_pred_ub_50")


### Functions to approximately evaluate the predictive distribution from a surveillance fit from hhh4
### These functions are based closely on simulate.hhh4 and simHHH4 from the surveillance package, but
### we have replaced calls to functions that sample from a poisson or negative binomial distribution
### with calls to functions that evaluate the density of a poisson or negative binomial distribution.
### Basically, we simulate nsim trajectories of length h - 1, where h is the largest horizon for which
### we want to evaluate the predictive distribution.  The predictive distribution is approximated as
### a mixture of the predictive distributions from these nsim trajectories.

get_predictive_density_and_sample_hhh4 <- function (object, # result from a call to hhh4
                           nsim=1, # number of replicates to simulate
                           seed=NULL,
                           y.start=NULL, # initial counts for epidemic components
                           subset=1:nrow(object$stsObj),
                           coefs=coef(object), # coefficients used for simulation
                           components=c("ar","ne","end"), # which comp to include
                           simplify=nsim>1, # counts array only (no full sts)
                           ...)
{
  ## Determine seed (this part is copied from stats:::simulate.lm with
  ## Copyright (C) 1995-2012 The R Core Team)
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                     # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ## END seed
  
  cl <- match.call()
  theta <- if (missing(coefs)) coefs else checkCoefs(object, coefs)
  
  ## lags
  lag.ar <- object$control$ar$lag
  lag.ne <- object$control$ne$lag
  maxlag <- max(lag.ar, lag.ne)
  
  ## initial counts
  nUnits <- object$nUnit
  if (is.null(y.start)) { # set starting value to mean observed (in subset!)
    y.means <- ceiling(colMeans(observed(object$stsObj)[subset,,drop=FALSE]))
    y.start <- matrix(y.means, maxlag, nUnits, byrow=TRUE)
  } else {
    if (is.vector(y.start)) y.start <- t(y.start)
    if (ncol(y.start) != nUnits)
      stop(sQuote("y.start"), " must have nUnits=", nUnits, " columns")
    if (nrow(y.start) < maxlag)
      stop("need 'y.start' values for lag=", maxlag, " initial time points")
  }
  
  ## get fitted components nu_it (with offset), phi_it, lambda_it, t in subset
  model <- surveillance:::terms.hhh4(object)
  means <- surveillance:::meanHHH(theta, model, subset=subset)
  psi <- surveillance:::splitParams(theta,model)$overdisp
  
  ## weight matrix/array of the ne component
  neweights <- surveillance:::getNEweights(object, coefW(theta))
  
  ## set predictor to zero if not included ('components' argument)
  stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
  getComp <- function (comp) {
    sel <- if (comp == "end") "endemic" else paste(comp, "exppred", sep=".")
    res <- means[[sel]]
    if (!comp %in% components) res[] <- 0
    res
  }
  ar <- getComp("ar")
  ne <- getComp("ne")
  end <- getComp("end")
  
  ## simulate and evaluate distribution conditional on one trajectory
  simcall <- quote(
    dHHH4(ar, ne, end, psi, neweights, y.start, lag.ar, lag.ne, observed(object$stsObj)[subset])
  )
  if (!simplify) {
    ## result template
    res0 <- object$stsObj[subset,]
    setObserved <- function (observed) {
      res0@observed[] <- observed
      res0
    }
    simcall <- call("setObserved", simcall)
  }
  # res <- if (nsim==1) eval(simcall) else
  #   replicate(nsim, eval(simcall),
  #             simplify=if (simplify) "array" else FALSE)
  temp <- lapply(seq_len(nsim), function(rep_ind) {
    eval(simcall)
  })
  
  res <- list(
    y_obs = object$stsObj[subset,],
    y_sampled = sapply(
      seq_len(nsim),
      function(rep_ind) {
        temp[[rep_ind]]$y_sampled
      },
      simplify = "array"
    ),
    obs_log_dens = sapply(
      seq_len(nsim),
      function(rep_ind) {
        temp[[rep_ind]]$obs_log_dens
      },
      simplify = "array"
    )
  )
  
  if(is.null(dim(res$y_sampled))) {
    dim(res$y_sampled) <- c(1, length(res$y_sampled))
  }
  if(is.null(dim(res$obs_log_dens))) {
    dim(res$obs_log_dens) <- c(1, 1, length(res$obs_log_dens))
  }
  
  # if (simplify) {
  #   dimnames(res)[1:2] <- list(subset, colnames(model$response))
  #   attr(res, "initial") <- y.start
  #   attr(res, "stsObserved") <- object$stsObj[subset,]
  #   class(res) <- "hhh4sims"
  # }
  
  ## Done
  attr(res, "call") <- cl
  attr(res, "seed") <- RNGstate
  res
}


dHHH4 <- function(ar,     # lambda_it (nTime x nUnits matrix)
                    ne,     # phi_it (nTime x nUnits matrix)
                    end,    # nu_it (nTime x nUnits matrix, offset included)
                    psi,    # overdisp param(s) or numeric(0) (psi->0 = Poisson)
                    neW,    # weight matrix/array for neighbourhood component
                    start,  # starting counts (vector of length nUnits, or
                    # matrix with nUnits columns if lag > 1)
                    lag.ar = 1,
                    lag.ne = lag.ar,
                  y_obs
)
{
  nTime <- nrow(end)
  nUnits <- ncol(end)
  
  ## simulate from Poisson or NegBin model
  rdistr <- if (length(psi)==0 ||
                isTRUE(all.equal(psi, 0, check.attributes=FALSE))) {
    rpois
  } else {
    psi.inv <- 1/psi   # since R uses different parametrization
    ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
    ## and overdispersion psi such that Variance = mean + psi*mean^2
    ## where 'size'=1/psi and length(psi) == 1 or length(mean)
    function(n, mean) rnbinom(n, mu = mean, size = psi.inv)
  }
  
  ## evaluate Poisson or NegBinmodel
  ddistr <- if (length(psi)==0 ||
                isTRUE(all.equal(psi, 0, check.attributes=FALSE))) {
    dpois
  } else {
    psi.inv <- 1/psi   # since R uses different parametrization
    ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
    ## and overdispersion psi such that Variance = mean + psi*mean^2
    ## where 'size'=1/psi and length(psi) == 1 or length(mean)
    function(x, mean, log) dnbinom(x, mu = mean, size = psi.inv, log = log)
  }
  
  ## if only endemic component -> simulate independently
  if (all(ar + ne == 0)) {
    return(matrix(rdistr(length(end), end), nTime, nUnits))
  }
  
  ## weighted sum of counts of other (neighbouring) regions
  ## params: y - vector with (lagged) counts of regions
  ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
  wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
    function (y, W) numeric(nUnits)
  } else function (y, W) .colSums(W * y, nUnits, nUnits)
  
  ## initialize matrices for means mu_i,t simulated data y_i,t and evaluated predictive density at observed values
  mu <- y <- obs_log_dens <- matrix(0, nTime, nUnits)
  y <- rbind(start, y)
  nStart <- nrow(y) - nrow(mu)        # usually just 1 for lag=1
  
  ## simulate and evaluate predictive density
  timeDependentWeights <- length(dim(neW)) == 3
  if (!timeDependentWeights) neWt <- neW
  for(t in seq_len(nTime)){
    if (timeDependentWeights) neWt <- neW[,,t]
    ## mean mu_i,t = lambda*y_i,t-1 + phi*sum_j wji*y_j,t-1 + nu_i,t
    mu[t,] <-
      ar[t,] * y[nStart+t-lag.ar,] +
      ne[t,] * wSumNE(y[nStart+t-lag.ne,], neWt) +
      end[t,]
    ## Sample from Poisson/NegBin with that mean
    y[nStart+t,] <- rdistr(nUnits, mu[t,])
    ## evaluate Poisson/NegBin with that mean at the observed value
    obs_log_dens[t,] <- ddistr(y_obs[t], mu[t, ], log = TRUE)
  }
  
#  ## return simulated data without initial counts
#  y[-seq_len(nStart),,drop=FALSE]
  ## return log predictive density of observed values
  return(list(y_sampled = y[- seq_len(nStart), drop = FALSE], obs_log_dens = obs_log_dens))
}



### Load data set and set variables describing how the fit is performed

## Load data for Dengue fever in San Juan
data <- read.csv("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/data-raw/San_Juan_Testing_Data.csv")

## Row indices in data corresponding to times at which we want to make a prediction
prediction_time_inds <- which(data$season %in% paste0(2009:2012, "/", 2010:2013))

## convert dates
data$time <- ymd(data$week_start_date)

## Add time_index column.  This is used for calculating the periodic kernel.
## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
## The origin is arbitrary.
data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))


## load surveillance fits and choose best one
surveillance_fits <- readRDS(file = file.path(
  "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
  data_set,
  "estimation-results/surveillance-fits.rds"))

aic_by_surveillance_fit <- sapply(surveillance_fits$model_fits, function(sfit) {
  summary(sfit)$AIC
})
bic_by_surveillance_fit <- sapply(surveillance_fits$model_fits, function(sfit) {
  summary(sfit)$BIC
})

if(identical(selection_criteria, "log_score")) {
  best_spec_ind <- which.min(surveillance_fits$model_specifications$mean_log_score)
  surveillance_fit <- surveillance_fits$model_fits[[best_spec_ind]]
} else if(identical(selection_criteria, "aic")) {
  best_spec_ind <- which.min(aic_by_surveillance_fit)
  surveillance_fit <- surveillance_fits$model_fits[[best_spec_ind]]
} else if(identical(selection_criteria, "bic")) {
  best_spec_ind <- which.min(bic_by_surveillance_fit)
  surveillance_fit <- surveillance_fits$model_fits[[best_spec_ind]]
}

n_sims <- 10000

# num_rows <- length(all_prediction_horizons) *
#   length(prediction_time_inds)
num_rows <- 1L

data_set_results <- data.frame(data_set = data_set,
                               model = "SARIMA",
                               prediction_horizon = rep(NA_integer_, num_rows),
                               prediction_time = rep(NA, num_rows),
                               log_score = rep(NA_real_, num_rows),
                               pt_pred = rep(NA_real_, num_rows),
                               AE = rep(NA_real_, num_rows),
                               interval_pred_lb_95 = rep(NA_real_, num_rows),
                               interval_pred_ub_95 = rep(NA_real_, num_rows),
                               interval_pred_lb_50 = rep(NA_real_, num_rows),
                               interval_pred_ub_50 = rep(NA_real_, num_rows),
                               stringsAsFactors = FALSE
)
class(data_set_results$prediction_time) <- class(data$time)

data_set_results <- data_set_results[-1, ]

# results_row_ind <- 1L

all_analysis_time_inds <- NULL
for(prediction_horizon in all_prediction_horizons) {
  for(prediction_time_ind in prediction_time_inds) {
    all_analysis_time_inds <- c(all_analysis_time_inds, prediction_time_ind - prediction_horizon)
  }
}

all_analysis_time_inds <- sort(unique(all_analysis_time_inds))

for(analysis_time_ind in all_analysis_time_inds) {
  ## Set values describing case in data_set_results
  # data_set_results$prediction_horizon[results_row_ind] <-
  #   prediction_horizon
  # data_set_results$prediction_time[results_row_ind] <-
  #   data$time[prediction_time_ind]
  
  max_prediction_horizon <- min(52L, nrow(data) - analysis_time_ind)
  
  
  ## Compute log score of distribution prediction at each horizon up to max_prediction_horizon
  predictive_density_and_sample <- get_predictive_density_and_sample_hhh4(surveillance_fit,
    nsim = n_sims,
    y.start = data$total_cases[analysis_time_ind],
    subset = seq(from = analysis_time_ind + 1, to = analysis_time_ind + max_prediction_horizon),
    simplify = TRUE)
  log_scores_by_sim_and_time <- predictive_density_and_sample$obs_log_dens[, 1, ] - log(n_sims)
  if(is.null(dim(log_scores_by_sim_and_time))) {
    dim(log_scores_by_sim_and_time) <- c(1, length(log_scores_by_sim_and_time))
  }
  log_scores_by_prediction_time <- logspace_sum_matrix_rows(log_scores_by_sim_and_time)
  
  phs_in_prediction_times <- which((analysis_time_ind + 1:52) %in% prediction_time_inds)
  prediction_time_inds_to_keep <- analysis_time_ind + phs_in_prediction_times
  
  new_results <- data.frame(data_set = data_set,
                            model = "HHH4",
                            prediction_horizon = phs_in_prediction_times,
                            prediction_time = data$time[prediction_time_inds_to_keep],
                            log_score = log_scores_by_prediction_time[phs_in_prediction_times],
                            pt_pred = rep(NA_real_, num_rows),
                            AE = rep(NA_real_, num_rows),
                            interval_pred_lb_95 = rep(NA_real_, num_rows),
                            interval_pred_ub_95 = rep(NA_real_, num_rows),
                            interval_pred_lb_50 = rep(NA_real_, num_rows),
                            interval_pred_ub_50 = rep(NA_real_, num_rows),
                            stringsAsFactors = FALSE
  )
#  class(data_set_results$prediction_time) <- class(data$time)
  
  y_sampled <- predictive_density_and_sample$y_sampled
  if(is.null(dim(y_sampled))) {
    dim(y_sampled) <- c(1, length(y_sampled))
  }
  
  ## Compute point prediction
  new_results$pt_pred <- apply(y_sampled[phs_in_prediction_times, , drop = FALSE], 1, median)

  ## Compute absolute error of point prediction
  new_results$AE <-
    abs(new_results$pt_pred -
      predictive_density_and_sample$y_obs@observed[phs_in_prediction_times])
  
  ## Compute prediction interval bounds
  new_results[, c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
    t(apply(y_sampled[phs_in_prediction_times, , drop = FALSE], 1, quantile, probs = c(0.025, 0.975, 0.25, 0.75)))
  
  ## append new_results to data_set_results
  data_set_results <- rbind(data_set_results, new_results)
} # analysis_time_ind

## Save results for the given data set
saveRDS(data_set_results, file = file.path(
  "/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
  data_set,
  "prediction-results/surveillance-predictions.rds"))
