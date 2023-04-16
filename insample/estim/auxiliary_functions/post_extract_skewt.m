function [samples_mu_coef, samples_sigma_coef, samples_alpha_coef, samples_nu, samples_mu, samples_sigma, samples_alpha] = post_extract_skewt(post_samples, warmup, thin, W, X)

    % Extract draws of skew-t parameter coefficients and time paths

    % Coefficients
    samples_mu_coef = post_extract_gambet(post_samples, 'mu', warmup, thin);
    samples_sigma_coef = post_extract_gambet(post_samples, 'sigma', warmup, thin);

    if isfield(post_samples(1), 'nu')
        samples_alpha_coef = post_extract_gambet(post_samples, 'alpha', warmup, thin);
        samples_nu = post_extract(post_samples, 'nu', warmup, thin);
    else
        samples_alpha_coef = zeros(size(samples_mu_coef));
        samples_nu = Inf(size(samples_mu_coef,1),1);
    end
    
    % Time paths
    samples_mu = samples_mu_coef*[W X]';
    samples_sigma = exp(samples_sigma_coef*[W X]');
    samples_alpha = samples_alpha_coef*[W X]';

end

function samples = post_extract_gambet(all_samples_struct, param, warmup, thin)

    % Extract samples of gamma and beta from all chains, after warm-up
    samples_gamma = post_extract(all_samples_struct, strcat('gamma_', param), warmup, thin);
    samples_beta = post_extract(all_samples_struct, strcat('beta_', param), warmup, thin);
    samples = [samples_gamma samples_beta];

end