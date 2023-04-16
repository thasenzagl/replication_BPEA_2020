% Estimate dynamic model using Stan
% MPM 2020-01-19

% Requires "settings" struct and "dat" table

addpath('auxiliary_functions/dfm');


%% Log

mkdir(settings.output_folder); % Create output folder if it doesn't already exist
delete(fullfile(settings.output_folder, 'mcmc.log'));
diary(fullfile(settings.output_folder, 'mcmc.log')); % Retain all screen output in log file


%% Data

% Read variable names and dates
varnames = dat.Properties.VariableNames;
date = datenum(dat.(settings.date_var));

% LHS variable Y
Y = dat.(settings.y_var);

% Unpenalized covariates
W = [[Y(1:end-1); nan] ones(length(Y),1)]; % Unpenalized covariates: lagged Y and const
w_vars = {'ylag', 'const'};

% Select penalized covariates
x_vars_ind = ~strcmp(varnames,settings.y_var) & ~strcmp(varnames,settings.date_var);
x_vars = varnames(x_vars_ind);
X = dat{:,x_vars};

% Impute missing observations and standardize X
if ~isfield(settings, 'r_imp') || ~isfield(settings, 'tol_imp')
    settings.r_imp = [];
    settings.tol_imp = [];
end
X = dfm_impute(X, settings.r_imp, settings.tol_imp);

% Select final sample, taking lag into account
Y = Y(2:end);
W = W(1:end-1,:);
X = X(1:end-1,:);

% Dimensions
T = length(Y);
q = size(W,2);
p = size(X,2);


%% Preliminary coefficient estimates

Xt = [W X];
estim_mu = Xt\Y; % Location
res = Y-Xt*estim_mu;
estim_sigma = Xt\(log(abs(res))+0.64); % Scale; note that E[log|N(0,1)|]=-0.64
estim_mu_scale = std(estim_mu(q+1:end));
estim_sigma_scale = std(estim_sigma(q+1:end));


%% Run Stan

% Initial parameter values
stan_init = struct('gamma_mu', estim_mu(1:q), ...
                   'beta_mu_std', estim_mu(q+1:end)/estim_mu_scale, ...
                   'gamma_sigma', estim_sigma(1:q), ...
                   'beta_sigma_std', estim_sigma(q+1:end)/estim_sigma_scale, ...
                   'lambda_mu', ones(p,1), ...
                   'lambda_sigma', ones(p,1), ...
                   'tau_mu', estim_mu_scale, ...
                   'tau_sigma', estim_sigma_scale);

if strcmp(settings.model, 'skewt')
    stan_init.gamma_alpha = 0.01*randn(q,1);
    stan_init.beta_alpha = 0.01*randn(p,1);
    stan_init.nu = 10;
end

% Data structure
stan_data = struct('T', T, ...
                   'q', q, ...
                   'p', p, ...
                   'Y', Y, ...
                   'W', {{W}}, ...
                   'X', {{X}});

if strcmp(settings.model, 'condhet')
    stan_data.min_scale = 0.1/p;
end

% Fit Stan model
model_prior = settings.model;
if isfield(settings, 'prior') && strcmp(settings.prior, 'horseshoe')
    model_prior = strcat(model_prior, '_horseshoe');
end
stan_fit = stan('file', strcat(model_prior, '.stan'), ...
                'data', stan_data, ...
                'init', stan_init, ...    
                'working_dir', settings.output_folder, ...
                settings.stan);
stan_fit.block();

% Print results
print(stan_fit);

% Save
post_samples = stan_fit.extract('permuted', false, 'inc_warmup', true); % Extract posterior samples (including warm-up)
save(fullfile(settings.output_folder, 'results.mat'), ...
     'settings', 'dat', 'date', 'w_vars', 'x_vars', 'Y', 'X', 'W', 'T', 'q', 'p', 'estim_*', ...
     'stan_init', 'stan_data', 'post_samples');


diary off;