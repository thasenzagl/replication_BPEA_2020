data {
    int<lower=1> T;         // No. of time periods
    int<lower=0> q;         // No. of unpenalized regressors
    int<lower=0> p;         // No. of penalized regressors
    vector[T] Y;            // Outcomes
    matrix[T,q] W;          // Unpenalized regressors
    matrix[T,p] X;          // Penalized regressors
    real<lower=0> min_scale;// Minimum for scale parameters;
}

parameters {
    vector[q] gamma_mu;             // Unpenalized coefficients: location
    vector[q] gamma_sigma;          // Unpenalized coefficients: scale
    vector[p] beta_mu_std;          // Penalized coefficients: location (standardized)
    vector[p] beta_sigma_std;       // Penalized coefficients: scale (standardized)
    vector[p] lambda_mu;            // Individual scale hyperparameters: location
    vector[p] lambda_sigma;         // Individual scale hyperparameters: scale
    real<lower=min_scale> tau_mu;   // Overall scale hyperparameter: location
    real<lower=min_scale> tau_sigma;// Overall scale hyperparameter: scale
}

transformed parameters {
    vector[p] beta_mu = tau_mu * (lambda_mu .* beta_mu_std);
    vector[p] beta_sigma = tau_sigma * (lambda_sigma .* beta_sigma_std);
}

model {
    // Prior
    gamma_mu ~ cauchy(0, 5);
    gamma_sigma ~ cauchy(0, 5);
    beta_mu_std ~ normal(0, 1);
    beta_sigma_std ~ normal(0, 1);
    lambda_mu ~ cauchy(0, 1);
    lambda_sigma ~ cauchy(0, 1);
    tau_mu ~ cauchy(0, 1);
    tau_sigma ~ cauchy(0, 1);

    // Likelihood
    Y ~ normal(W*gamma_mu + X*beta_mu, exp(W*gamma_sigma + X*beta_sigma));
}

