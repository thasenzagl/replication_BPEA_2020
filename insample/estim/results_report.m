%clear;
addpath('auxiliary_functions/', 'auxiliary_functions/skew_t/');

% Produce figures and tables of MCMC results
% MPM 2020-01-20


%% Settings

% Specification
%spec = 'us';        % Either 'us' or 'global'
%predic = 'factor';  % Either 'factor' or 'var'
%start_year = [];    % Start year (set to empty [] if full sample)

% Model
%model = 'skewt';    % Either 'condhet' or 'skewt'

% Files
results_root = 'output/stan';   % Root folder with MCMC results
save_root = 'output/fig_tab';   % Root folder for saving figures/tables
save_eps = true;                % true: save density plots as .eps files (in addition to .png)

% Plot/table settings
%zero_factor = [];                       % Index j of x_j to set equal to 0 when plotting time-varying parameters/moments (set to [] if actual data should be used)
plot_xlim_beta_sigma = [-0.5 0.5];      % Horiz axis limits for beta_sigma plots
plot_param_quants = [0.25 0.75];        % Quantiles to mark in parameter plots
plot_moment_quants = [0.05 0.5 0.95];   % Three quantiles to report in time-varying moment plots (second one is marked with thick line)
tab_quants = plot_param_quants;         % Quantiles to store in table
magn_thresh_mu = 0.2/2;                 % Report posterior probability of mu coefs exceeding this threshold in abs. val.
magn_thresh_sigma = 0.1/2;              % Report posterior probability of sigma coefs exceeding this threshold in abs. val.
magn_thresh_tvdape = 0.05/2;            % Report posterior probability of TVD APEs exceeding this threshold in abs. val.


%% Find results folders

results_folder = fullfile(spec, predic, strcat(model, num2str(start_year))); % EITHER: (i) subfolder with results.mat file, OR (ii) subfolder with several sub-subfolders that contain results.mat files
results_subfolder = fullfile(results_root, results_folder); % Results subfolder
save_subfolder = fullfile(save_root, results_folder); % Save subfolder

if isfile(fullfile(results_subfolder, 'results.mat')) % If results.mat file exists here...
    results_subsubfolders = {''}; % Just this folder
else % Otherwise, scan through subfolders...
    the_dir = dir(results_subfolder); % All contents of this folder
    results_subsubfolders = {the_dir([false false [the_dir(3:end).isdir]==1]).name}; % All subfolders of this folder (except '.' and '..')    
end

% Types of parameters to summarize
param_types = {'mu','sigma'};
if strcmp(model, 'skewt')
    param_types = [param_types {'alpha'}];
end


numsubsubfolders = length(results_subsubfolders);
results_table = table;

for i=1:numsubsubfolders % For each results subfolder...
    
    disp(results_subsubfolders{i});
    
    % Save subsubfolder
    the_save_subsubfolder = fullfile(save_subfolder, results_subsubfolders{i});
    mkdir(the_save_subsubfolder); % Create save folder
    
    
    %% Load posterior samples
    
    load(fullfile(results_subfolder, results_subsubfolders{i}, 'results.mat'));
    the_warmup = settings.stan.warmup;
    the_thin = settings.stan.thin;
    
    
    %% Figures: trace plots and posterior densities of model parameters

    for g={'gamma','beta'}
        
        % Determine plot font size
        if strcmp(g, 'gamma') || strcmp(predic, 'factor')
            the_fontsize = 12;
        else
            the_fontsize = 6;
        end

        for s=param_types

            the_param = strcat(g{1}, '_', s{1});
            
            % Extract parameter samples for each chain
            [the_samples_nowarmup, the_samples] = post_extract(post_samples, the_param, the_warmup, the_thin);

            np = size(the_samples,3); % Number of parameter components
            numrow = floor(sqrt(np)); % Number of rows in plot
            numcol = ceil(np/numrow); % Number of columns in plot
            f_trace=figure;
            f_dens=figure;

            for j=1:np % For each component in parameter vector...
                
                % Title of subplot
                if strcmp(g{1}, 'gamma')
                    the_title = w_vars{j};
                else
                    the_title = x_vars{j};
                end

                % TRACE PLOT
                
                figure(f_trace);
                subplot(numrow,numcol,j);
                
                plot((1:size(the_samples,2))*the_thin, the_samples(:,:,j), '-'); % One line for each chain
                the_ylim = ylim;
                hold on;
                line(the_warmup*ones(1,2), the_ylim, 'Color', 'k', 'LineStyle', '--'); % Add vertical line showing warm-up period finish
                hold off;
                ylim(the_ylim);
                set(gca, 'FontSize', the_fontsize);
                title(the_title, 'FontSize', the_fontsize, 'Interpreter', 'none');
                
                % POSTERIOR DENSITY PLOT
                
                figure(f_dens);
                subplot(numrow,numcol,j);
                
                [the_fval,the_xi] = ksdensity(the_samples_nowarmup(:,j)); % Smoothed posterior density
                the_quants = quantile(the_samples_nowarmup(:,j), plot_param_quants); % Quantiles

                plot(the_xi, the_fval, '-k', 'LineWidth', 1);

                % Add vertical axis and quantiles
                the_ylim = ylim;
                hold on;
                plot([0 0], the_ylim, 'Color', 'k', 'LineStyle', '-');
                plot(the_quants(1)*ones(1,2), the_ylim, 'Color', 'r', 'LineStyle', '--');
                plot(the_quants(2)*ones(1,2), the_ylim, 'Color', 'r', 'LineStyle', '--');
                hold off;
                ylim(the_ylim);

                % Horizontal axis
                if strcmp(g{1}, 'beta') && strcmp(s{1}, 'sigma')
                    xlim(plot_xlim_beta_sigma);
                end

                set(gca, 'FontSize', the_fontsize);
                title(the_title, 'FontSize', the_fontsize, 'Interpreter', 'none');

            end

            % Save figures
            save_fig(f_trace, fullfile(the_save_subsubfolder, strcat('trace_', the_param)), save_eps);
            save_fig(f_dens, fullfile(the_save_subsubfolder, strcat('postdens_', the_param)), save_eps);

        end

    end


    %% Table: summary statistics
    
    the_X_plot = X;
    if ~isempty(zero_factor)
        the_X_plot(:,zero_factor)=0; % Set x_j equal to 0 at all points in time
    end
    
    % Extract skew-t parameter coefficients and time paths
    [samples_mu_coef, samples_sigma_coef, samples_alpha_coef, samples_nu, samples_mu, samples_sigma, samples_alpha] = post_extract_skewt(post_samples, the_warmup, the_thin, W, the_X_plot);
    [numdraw,T] = size(samples_mu);
    
    % Compute posterior summaries for mu and sigma
    [probpos_mu, probverypos_mu, probveryneg_mu, med_mu, qlo_mu, qhi_mu] = post_summ(samples_mu_coef, magn_thresh_mu, tab_quants);
    [probpos_sigma, probverypos_sigma, probveryneg_sigma, med_sigma, qlo_sigma, qhi_sigma] = post_summ(samples_sigma_coef, magn_thresh_sigma, tab_quants);
    
    % Gather in results table
    var = [w_vars x_vars]';
    numvars = length(var);
    geo = repmat(results_subsubfolders(i), numvars, 1);
    the_table = table(geo, var, ...
                      probpos_mu, probverypos_mu, probveryneg_mu, med_mu, qlo_mu, qhi_mu, ...
                      probpos_sigma, probverypos_sigma, probveryneg_sigma, med_sigma, qlo_sigma, qhi_sigma);
    
    if strcmp(model, 'skewt')
        
        % Compute average partial effects on TVD = atan(abs(alpha))/pi
        samples_tvd = atan(abs(samples_alpha))/pi; % Samples of TVD(alpha_t)
        samples_tvdavgderiv = mean(sign(samples_alpha)./(1+samples_alpha.^2),2)/pi; % Samples of avg_t{dTVD/dalpha_t}
        samples_tvdape = samples_tvdavgderiv.*samples_alpha_coef; % Samples of avg_t{dTVD/dx_{jt}}
        
        % Compute posterior summaries of TVD APEs
        [probpos_alpha, probverypos_alpha, probveryneg_alpha, med_alpha, qlo_alpha, qhi_alpha] = post_summ(samples_tvdape, magn_thresh_tvdape, tab_quants);
        
        % Compute posterior summaries of alpha and TVD
        mean_alphaidx = repmat(mean(samples_alpha(:)), numvars, 1);     % Mean time-average of alpha
        mean_tvdavg = repmat(mean(samples_tvd(:)), numvars, 1);         % Mean time-average of TVD
        mean_tvdstd = repmat(mean(std(samples_tvd,0,2)), numvars, 1);   % Mean time-stdev of TVD
        
        % Summarize nu (d.f. parameter)
        [~, ~, ~, med_nu, qlo_nu, qhi_nu] = post_summ(samples_nu, 0, tab_quants); % Summary stats
        med_nu = repmat(med_nu, numvars, 1);
        qlo_nu = repmat(qlo_nu, numvars, 1);
        qhi_nu = repmat(qhi_nu, numvars, 1);
        
        % Add to results table
        the_table_alpha = table(probpos_alpha, probverypos_alpha, probveryneg_alpha, med_alpha, qlo_alpha, qhi_alpha, mean_alphaidx, mean_tvdavg, mean_tvdstd, med_nu, qlo_nu, qhi_nu);
        the_table = [the_table the_table_alpha];
        
    end

    % Append table to existing results
    results_table = [results_table; the_table];
    
    
    %% Figures: conditional moments over time
    
    zero_factor_suffix = '';
    if ~isempty(zero_factor)
        zero_factor_suffix = sprintf('%s%d', '_nofactor', zero_factor);
    end
    
    % PLOT: PARAMETER EVOLUTION
    figure;
    plot_ts_quant_band(date(1:end-1)', {samples_mu, samples_sigma, samples_alpha, repmat(samples_nu,1,T)}, plot_moment_quants, {'mu', 'sigma', 'alpha', 'nu'});
    save_fig(gcf, fullfile(the_save_subsubfolder, strcat('ts_param', zero_factor_suffix)), save_eps);

    % Compute implied conditional moments
    [samples_mean, samples_var, samples_skew, samples_kurt] = skew_t_moments(samples_alpha, samples_nu); % Moments, assuming mu=0, sigma=1
    samples_mean = samples_mu + samples_sigma.*samples_mean; % Adjust for mu and sigma
    samples_std = samples_sigma.*sqrt(samples_var); % Adjust for sigma
    
    % Attempt to load pre-computed results
    the_compute_filename = fullfile(the_save_subsubfolder, 'compute.mat');
    if isempty(zero_factor) && isfile(the_compute_filename)
        load(the_compute_filename);
    else
        h = [];
        lower_quant = [];
    end
    
    for l=[1 h] % For each forecast horizon...

        % PLOT: CONDITIONAL MOMENT EVOLUTION
        figure;
        if l==1
            the_curves = {4*samples_mean, 4*samples_std, samples_skew, samples_kurt}; % Annualized mean and std dev
        else
            the_curves = {4/h*samples_mean_h, 4/h*samples_std_h, samples_skew_h, samples_kurt_h};
        end
        plot_ts_quant_band(date(1:end-1)', the_curves, plot_moment_quants, {'mean', 'standard deviation', 'skewness', 'kurtosis'});
        save_fig(gcf, fullfile(the_save_subsubfolder, sprintf('%s%d%s', 'ts_moments_h', l, zero_factor_suffix)), save_eps);

        % PLOT: CONDITIONAL PROBABILITY/SHORTFALL EVOLUTION
        if ~isempty(lower_quant)
            figure;
            if l==1
                samples_shortfall_minusmean = samples_shortfall - samples_mean; % Expected shortfall minus mean
                the_curves = {samples_prob_neg, samples_prob_belowmean, 4*samples_shortfall, 4*samples_shortfall_minusmean};
            else
                samples_shortfall_minusmean_h = samples_shortfall_h - h*samples_mean(draws_h,:);
                the_curves = {samples_prob_neg_h, samples_prob_belowmean_h, 4/h*samples_shortfall_h, 4/h*samples_shortfall_minusmean_h};
            end
            plot_ts_quant_band(date(1:end-1)', the_curves, plot_moment_quants, {'probability of negative growth', 'probability of growth below mean', sprintf('%d%s', 100*lower_quant, '% expected shortfall'), sprintf('%d%s', 100*lower_quant, '% expected shortfall minus mean')});
            save_fig(gcf, fullfile(the_save_subsubfolder, sprintf('%s%d', 'ts_other_h', l)), save_eps);
        end

    end
    
end


%% Save table to CSV

writetable(results_table, fullfile(save_subfolder, 'results.csv'));


%% Auxiliary functions for posterior summaries

function [probpos, probverypos, probveryneg, med, qlo, qhi] = post_summ(samples, magn_thresh, tab_quants)

    % Compute summary stats
    probpos = mean(samples > 0)';
    probverypos = mean(samples > magn_thresh)';
    probveryneg = mean(samples < -magn_thresh)';
    med = median(samples)';
    qlo = quantile(samples, min(tab_quants))';
    qhi = quantile(samples, max(tab_quants))';

end

function plot_ts_quant_band(dates, data, quants, titles)

    % Posterior quantiles as time series plot

    % Compute quantiles for each time-varying quantity
    numplot = length(data);
    curves = cell(1,length(numplot));
    for j=1:numplot
        curves{j} = quantile(data{j}, quants);
    end
    
    % Plot
    plot_ts_band(dates, curves, titles);

end

function plot_ts_band(dates, curves, titles)

    % Plot time series with shaded band

    numplot = length(curves);
    
    for j=1:numplot
        
        the_curves = curves{j};
        subplot(numplot,1,j);
        fill([dates, fliplr(dates)], [the_curves(1,:), fliplr(the_curves(3,:))], [0.7 0.7 0.7]); % Shaded band
        hold on;
        plot(dates, the_curves(2,:), '-k', 'LineWidth', 2); % Thick curve
        hold off;
        datetick('x');
        title(titles{j}, 'Interpreter', 'none');
        
    end

end

function save_fig(f, filename, save_eps)

    % Save figure as .png and possibly .eps

    saveas(f, strcat(filename, '.png'));
    if save_eps
        saveas(f, strcat(filename, '.eps'), 'epsc');
    end
    close(f);

end
