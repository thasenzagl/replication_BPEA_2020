data {
    int<lower=1> T;     // No. of time periods
    int<lower=0> q;     // No. of unpenalized regressors
    int<lower=0> p;     // No. of penalized regressors
    vector[T] Y;        // Outcomes
    matrix[T,q] W;      // Unpenalized regressors
    matrix[T,p] X;      // Penalized regressors
}

parameters {
    vector[q] gamma_mu;         // Unpenalized coefficients: location
    vector[q] gamma_sigma;      // Unpenalized coefficients: scale
    vector[p] beta_mu_std;      // Penalized coefficients: location (standardized)
    vector[p] beta_sigma_std;   // Penalized coefficients: scale (standardized)
    real<lower=0> tau_mu;       // Scale hyperparameter: location
    real<lower=0> tau_sigma;    // Scale hyperparameter: scale
}

transformed parameters {
    vector[p] beta_mu = tau_mu*beta_mu_std;
    vector[p] beta_sigma = tau_sigma*beta_sigma_std;
}

model {
    // Prior
    gamma_mu ~ cauchy(0, 5);
    gamma_sigma ~ cauchy(0, 5);
    beta_mu_std ~ normal(0, 1);
    beta_sigma_std ~ normal(0, 1);
    tau_mu ~ cauchy(0, 1);
    tau_sigma ~ cauchy(0, 1);

    // Likelihood
    Y ~ normal(W*gamma_mu + X*beta_mu, exp(W*gamma_sigma + X*beta_sigma));
}

