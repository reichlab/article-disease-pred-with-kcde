data {
    int<lower=1> N; // number of observations
    int<lower=1> M; // number of models being evaluated
    int<lower=1> B; // number of columns in spline basis for each model
    matrix[N, M * B] X; // matrix with spline basis (per model; M * B columns)
    real prediction_time[N]; // 
//    real prediction_horizon[N]; // 
//    real model[N]; // 
    real log_score_difference[N]; // outcome 
}
parameters {
    vector[M * B] beta; // coefficients for splines in prediction horizon for each model
    vector<lower=0>[M * B] sigma_h_sq; // prior variance for betas
    real<lower=0> sigma; // standard deviation parameter for errors around spline
}
transformed parameters {
    vector[N] mu; // means from splines
    mu <- X*beta;
}
model {
    // priors
    for(h in 1:(M * B)) {
        sigma_h_sq[h] ~ inv_gamma(0.25, 0.25);
        beta[h] ~ normal(0, sqrt(sigma_h_sq[h]));
    }
    sigma ~ inv_gamma(0.5, 0.5);
//    for(h in 1:(M * B)) {
//        beta[h] ~ cauchy(0, 2.5);
//    }
    
    // iid normal model
    log_score_difference ~ normal(mu, sigma);
}
generated quantities {
    vector[N] y_pred;
    y_pred <- X*beta; //the y values predicted by the model
}

