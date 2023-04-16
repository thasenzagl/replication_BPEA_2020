function [y_sim, x_sim] = sim_dyn_skewt_var(y, x, mu_coef, sigma_coef, alpha_coef, nu, A, Sigma)

    % Simulate one step from a dynamic skew-t model for y_t
    % and VAR(1) model for x_t
    
    % Simulate x
    x_sim = mvnrnd(A*x, Sigma)';
    
    % Compute skew-t parameters
    wx = [y; 1; x]; % Stack w_t and x_t
    mu = mu_coef*wx;
    sigma = exp(sigma_coef*wx);
    alpha = alpha_coef*wx;
    
    % Simulate y
    y_sim = mu + sigma*skew_t_sim(alpha, nu);

end