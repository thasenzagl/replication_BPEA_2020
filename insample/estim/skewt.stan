data {
    int<lower=1> T; // No. of time periods
    int<lower=0> q; // No. of unpenalized regressors
    int<lower=0> p; // No. of penalized regressors
    vector[T] Y;    // Outcomes
    matrix[T,q] W;  // Unpenalized regressors
    matrix[T,p] X;  // Penalized regressors
}

parameters {
    vector[q] gamma_mu;         // Unpenalized coefficients: location
    vector[q] gamma_sigma;      // Unpenalized coefficients: scale
    vector[q] gamma_alpha;      // Unpenalized coefficients: shape
    vector[p] beta_mu_std;      // Penalized coefficients: location (standardized)
    vector[p] beta_sigma_std;   // Penalized coefficients: scale (standardized)
    vector[p] beta_alpha_std;   // Penalized coefficients: shape (standardized)
    real<lower=0> nu;           // Degrees of freedom
    real<lower=0> tau_mu;       // Scale hyperparameter: location
    real<lower=0> tau_sigma;    // Scale hyperparameter: scale
	real<lower=0> tau_alpha;    // Scale hyperparameter: shape
}

transformed parameters {
    vector[p] beta_mu = tau_mu*beta_mu_std;
    vector[p] beta_sigma = tau_sigma*beta_sigma_std;
	vector[p] beta_alpha = tau_alpha*beta_alpha_std;
}

model {
    // Define auxiliary quantities for likelihood
    vector[T] mus;
    vector[T] log_sigmas;
    vector[T] alphas;
    vector[T] Y_std;
    vector[T] cdf_arg;

    // Prior
    gamma_mu ~ cauchy(0, 5);
    gamma_sigma ~ cauchy(0, 5);
    gamma_alpha ~ cauchy(0, 5);
    beta_mu_std ~ normal(0, 1);
    beta_sigma_std ~ normal(0, 1);
    beta_alpha_std ~ normal(0, 1);
    nu ~ gamma(1.5, 0.1);
    tau_mu ~ cauchy(0, 1);
    tau_sigma ~ cauchy(0, 1);
	tau_alpha ~ cauchy(0, 1);

    // Compute auxiliary quantities for likelihood
    mus = W*gamma_mu + X*beta_mu;
    log_sigmas = W*gamma_sigma + X*beta_sigma;
    alphas = W*gamma_alpha + X*beta_alpha;

    Y_std = (Y-mus).*exp(-log_sigmas);
    cdf_arg = sqrt((nu + square(Y_std)) / (nu + 1));

    // Increment with skew-t density
    target += - sum(log_sigmas)
              + student_t_lpdf(Y_std | nu, 0, 1) 
              + student_t_lcdf((alphas .* Y_std) ./ cdf_arg | nu+1, 0, 1);
}

