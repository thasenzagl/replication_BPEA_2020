clear;
addpath('auxiliary_functions/', 'auxiliary_functions/skew_t/', 'auxiliary_functions/var/');

% Computationally intensive posterior summaries of MCMC results
% MPM 2020-02-20

% ONLY WORKS FOR U.S. FACTOR SPECIFICATION!


%% Settings

% Specification
start_year = [];    % Start year (set to empty [] if full sample)

% Model
model = 'skewt';    % Either 'condhet' or 'skewt'

% Saved MCMC results
results_root = 'output/stan';       % Root folder with MCMC results
save_root = 'output/fig_tab';       % Root folder to save computations in

% Expected shortfall
lower_quant = 0.05;                 % Quantile used for computing expected shortfall

% h-step-ahead forecast
h = 4;                              % Forecast horizon (set to [] if h-step-ahead forecast not desired)
numdraw_h = 2e3;                    % Number of posterior parameter draws to consider
numsim_h = 5e3;                     % Number of forecast simulations per posterior parameter draw
rng(202002202);                     % Random number seed for forecast simulations

% Parallel computing
poolobj = parpool;                  % Open parallel pool


%% Load MCMC results

file_folder = fullfile('us', 'factor', strcat(model, num2str(start_year))); % Sub-folder with model/specification of interest
load(fullfile(results_root, file_folder, 'results.mat'));
warmup = settings.stan.warmup;
thin = settings.stan.thin;


%% Extract relevant posterior draws

% Extract skew-t parameter coefficients and time paths
[samples_mu_coef, samples_sigma_coef, samples_alpha_coef, samples_nu, samples_mu, samples_sigma, samples_alpha] = post_extract_skewt(post_samples, warmup, thin, W, X);

samples_mean = skew_t_moments(samples_alpha, samples_nu); % Conditional mean when mu=0, sigma=1
samples_mean = samples_mu + samples_sigma.*samples_mean; % Adjust for mu and sigma

[numdraw,T] = size(samples_mu);


%% h=1: Numerical computations

samples_prob_neg = nan(numdraw,T);
samples_prob_belowmean = nan(numdraw,T);
samples_shortfall = nan(numdraw,T);

opt_fsolve = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true); % Options for fsolve (when computing quantiles)

disp('h = 1: Computing probabilities, quantiles, shortfall numerically');
timer = tic;

parfor ip=1:numdraw % Loop over posterior draws

    for it=1:T % Loop over time

        % h=1 calculations
        the_cdf = skew_t_cdf(([0 samples_mean(ip,it)]-samples_mu(ip,it))/samples_sigma(ip,it), samples_alpha(ip,it), samples_nu(ip));
        samples_prob_neg(ip,it) = the_cdf(1); % Probability of negative outcome
        samples_prob_belowmean(ip,it) = the_cdf(2); % Probability of outcome below conditional mean
        the_quant = skew_t_quant(lower_quant, samples_alpha(ip,it), samples_nu(ip), opt_fsolve); % Quantile, assuming mu=0, sigma=1
        samples_shortfall(ip,it) = samples_mu(ip,it) + samples_sigma(ip,it)*skew_t_condmean(the_quant, samples_alpha(ip,it), samples_nu(ip)); % Expected shortfall

    end

    if mod(ip,ceil(numdraw/100))==0
        fprintf('%s%3d%s\n', repmat(' ',1,floor(50*ip/numdraw)), 100*ip/numdraw, '%');
    end

end

disp('Done. Elapsed time (min):');
disp(toc(timer)/60);


%% h steps ahead: Simulations

if ~isempty(h)

    [var_Ahat, var_Sigmahat, var_XpX, var_T] = var_estim(X, 1); % Estimate VAR for x
    var_Sigmahat_inv_chol = chol(inv(var_Sigmahat));

    draws_h = datasample(1:numdraw, numdraw_h, 'Replace', false); % Randomly select subset of posterior parameter draws

    samples_mu_coef_h = samples_mu_coef(draws_h,:);
    samples_sigma_coef_h = samples_sigma_coef(draws_h,:);
    samples_alpha_coef_h = samples_alpha_coef(draws_h,:);
    samples_nu_h = samples_nu(draws_h);
    
    samples_mean_h = nan(numdraw_h,T);
    samples_std_h = nan(numdraw_h,T);
    samples_skew_h = nan(numdraw_h,T);
    samples_kurt_h = nan(numdraw_h,T);
    samples_prob_neg_h = nan(numdraw_h,T);
    samples_prob_belowmean_h = nan(numdraw_h,T);
    samples_shortfall_h = nan(numdraw_h,T);

    fprintf('%s%d%s\n', 'h = ', h, ': Computing probabilities, quantiles, shortfall by simulation');
    timer = tic;

    parfor ip=1:numdraw_h % Loop over a smaller number of posterior draws

        % Draw from normal-inverse-Wishart posterior for VAR parameters
        [the_A_draw, the_Sigma_draw] = sim_post_niw_diffuse(var_Ahat, var_Sigmahat_inv_chol, var_XpX, var_T);

        for it=1:T % Loop over time

            % Simulate h steps
            sim_y_cum = zeros(numsim_h,1);
            for is=1:numsim_h
                the_y = W(it,1);
                the_x = X(it,:)';
                the_y_cum = 0;
                for l=1:h % Simulate forward one step at a time
                    [the_y, the_x] = sim_dyn_skewt_var(the_y, the_x, samples_mu_coef_h(ip,:), samples_sigma_coef_h(ip,:), samples_alpha_coef_h(ip,:), samples_nu_h(ip), the_A_draw, the_Sigma_draw);
                    the_y_cum = the_y_cum + the_y;
                end
                sim_y_cum(is) = the_y_cum; % Store simulated y_{t+1}+...+y_{t+h}
            end

            % Compute moments of y_{t+1}+...+y_{t+h}
            samples_mean_h(ip,it) = mean(sim_y_cum);
            samples_std_h(ip,it) = std(sim_y_cum);
            samples_skew_h(ip,it) = skewness(sim_y_cum);
            samples_kurt_h(ip,it) = kurtosis(sim_y_cum);

            % Compute other statistics
            samples_prob_neg_h(ip,it) = mean(sim_y_cum<0);
            samples_prob_belowmean_h(ip,it) = mean(sim_y_cum<h*samples_mean(ip,it));
            samples_shortfall_h(ip,it) = mean(sim_y_cum(sim_y_cum<quantile(sim_y_cum,lower_quant)));

        end

        if mod(ip,ceil(numdraw_h/100))==0
            fprintf('%s%3d%s\n', repmat(' ',1,floor(50*ip/numdraw_h)), 100*ip/numdraw_h, '%');
        end

    end

    disp('Done. Elapsed time (min):');
    disp(toc(timer)/60);

end

delete(poolobj); % Close parallel pool


%% Save results

% Variables to save
save_vars = {'lower_quant', 'h', 'samples_prob_neg', 'samples_prob_belowmean', 'samples_shortfall'};
if ~isempty(h)
    save_vars = [save_vars, {'numdraw_h', 'numsim_h', 'draws_h', 'var_*', 'samples_*_h'}];
end

save_subfolder = fullfile(save_root, file_folder);
mkdir(save_subfolder); % Create save folder
save(fullfile(save_subfolder, 'compute.mat'), save_vars{:}); % Save computations

