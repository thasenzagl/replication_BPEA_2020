clear all;
close all;
clc;

%%
% This script runs a real time forecasting exercise on real GDP
% growth using quantile regressions. First we use a real time dataset 
% and calendar from Alfred and construct a set of real time 
% vintages. Then, for each vintage we perform the following steps

%       1. We estimate monthly common factors and average the factors
%       across quarters
%       2. We estimate a quantile regression of the factors on real GDP
%       growth to obtain the predictive distributions for backcast,
%       nowcast, and forecast
%       3. We fit skewed student t distributions using the quantiles

% The script uses functions and code snippets that were written 
% for the following papers:
%       1. Banbura and Modugno (2010)
%       2. Adrien et al. (2019)

%% User Input

% Select specification (run script separately for different specifications)
% mode = 'global_and_fin': Run exercise with global and financial factors 
% mode = 'global_only': Run exercise with global factor only 
% mode = 'nonfin_only': Run exercise with nonfinancial factor only 
mode = 'global_and_fin';

% Sample Start date
P.StartEst = [1980 1];

% GDP Fred code
P.SerNews = 'GDPC1';

% Start and Enddates for pseudo real time vintage
StartDate = '01-Jan-2005';
EndDate   = '30-Sept-2019';

% Code Library
addpath(genpath('../../DFM_library/'));
addpath(genpath('../../QR_library/'));

% The number of lags in the VAR process of the factors. 
P.p=2;

% The maximum number of iterations in the EM algorithm
P.max_iter=1000;

% Annualize GDP Data
annualize = 1;

% Number of factors per block
P.nR = 1;

% Data
if strcmp(mode,'global_only') || strcmp(mode,'global_and_fin')
    P.Path = [pwd '/data/global_fin.xlsx'];
elseif strcmp(mode,'nonfin_only')
    P.Path = [pwd '/data/nonfin_only.xlsx'];
end


%% Quantile regression settings

deltaYY = 0.1;
YY = (-20):deltaYY:20;

QQ = 0.05:0.05:0.95;
% Find the location of the median in QQ
[~,jq50] = min(abs(QQ-.50));

% Target probability for expected shortfall/longrise
alpha = 0.05;

% Interval width for approximation to expected shortfall/longrise integral
delta1 = 0.01;

%% Initial settings for real and financial

% Inital Settings for estimation
[P, DatesV, DD, DD_series, name_vintage, StartVintage, EndVintage, blocks, target_last, include] = dfm_initial_settingsRT(P, StartDate, EndDate);

% Preallocate to store results
setup;

%% Run OOS forecasting exercise

for t = StartVintage:EndVintage

    disp("Now Running: "  + datestr(name_vintage(t)));
    
    % update real time data
    [P, iSer, iSer_idx]  = update_rt_data(P, DD_series{t}, blocks, include);
    
    %% Find iQ
    temp = datevec(name_vintage(t));
    temp = temp(1:2);
    
    % Turn months 1,2->3, 4,5->6, 7,8->9, 10,11->12
    if  mod(temp(2),3) == 1
        temp(2)= temp(2)+2;
    elseif mod(temp(2),3) == 2
        temp(2)= temp(2)+1;
    end
    
    P.Qnews = temp;
    
    % Index of DatesV where the date is the same as P.Qnews
    iQ = find(DatesV(:,1) == P.Qnews(1) & DatesV(:,2) == P.Qnews(2));
    
    %% Compute factor
    target = DD{t}(:,iSer);
    X_new_temp=DD{t}(:,iSer_idx); % remove GDP growth
    X_new = X_new_temp(:,P.include); % remove monthly variables that are not included
    
    if strcmp(mode,'global_only') || strcmp(mode,'global_and_fin')
        [fr, ff, y, datesQ, jt, first_not_na] = estim_factor(X_new, P, DatesV, iQ, target, annualize);
        factor_global(:,t) = fr;
        factor_fin(:,t) = ff; 
    elseif strcmp(mode,'nonfin_only')
        [fr, ~, y, datesQ, jt, first_not_na] = estim_factor(X_new, P, DatesV, iQ, target, annualize);
        factor_global(:,t) = fr;
    end    
        
    if strcmp(mode,'global_only') || strcmp(mode,'nonfin_only')
        Z = [ones(size(y)), factor_global(:,t)];
    elseif strcmp(mode,'global_and_fin')
        Z = [ones(size(y)), factor_global(:,t), factor_fin(:,t)];
    end    
        
    %% Compute forecasts for each quantile
    
    for jq = 1:length(QQ)
        
        % Backcast
        if isnan(y(jt-1))
            b = rq(Z(first_not_na:(jt - 2), :), y(first_not_na:(jt - 2)), QQ(jq));
            YQ_bc(t,jq) = Z(jt-1, :) * b;
            
            bunc = rq(Z(first_not_na:(jt - 2), 1), y(first_not_na:(jt - 2)), QQ(jq));
        else
            b = rq(Z(first_not_na:(jt - 1), :), y(first_not_na:(jt - 1)), QQ(jq));
            bunc = rq(Z(first_not_na:(jt - 1), 1), y(first_not_na:(jt - 1)), QQ(jq));
        end
        
        % Nowcast
        YQ_nc(t, jq) = Z(jt, :) * b;
        
        % Forecast
        YQ_fc(t, jq) = Z(jt+1, :) * b;
        
        % Unconditional
        YQunc(t,jq) = bunc;
        
    end
    
    %% Outturn
    if annualize == 1
        OT_bc(t) = 400 * target_last(iQ-3);
        OT_nc(t) = 400 * target_last(iQ);
        OT_fc(t) = 400 * target_last(iQ+3);
    else
        OT_bc(t) = target_last(iQ-3);
        OT_nc(t) = target_last(iQ);
        OT_fc(t) = target_last(iQ+3);
    end
    
    %% Fit skewed t-distribution
    
    % Unconditional
    [lc, sc, sh, df] = QuantilesInterpolation(YQunc(t, :), QQ);
    PSTunc(t,:) = dskt(YY, lc, sc, sh, df);
    
    % Backcast
    if isnan(y(jt-1))
        [lc, sc, sh, df] = QuantilesInterpolation(YQ_bc(t, :), QQ);
        Par_bc(t,:) = [lc, sc, sh, df];
        PST_bc(t,:) = dskt(YY, lc, sc, sh, df);
        QST_bc(t,:) = qskt(QQ, lc, sc, sh, df);
        PS_bc(t) = dskt(OT_bc(t), lc, sc, sh, df);
        
        qqTemp = qskt([delta1:delta1:alpha], Par_bc(t,1), Par_bc(t,2), Par_bc(t,3), Par_bc(t,4));
        SF_bc(t,:) = 1/alpha * sum(qqTemp * delta1);
    
        qqTemp = qskt([1-delta1:-delta1:1-alpha], Par_bc(t,1), Par_bc(t,2), Par_bc(t,3), Par_bc(t,4));
        LR_bc(t,:) = 1/alpha * sum(qqTemp * delta1);
    end
    
    % Nowcast
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_nc(t, :), QQ);
    Par_nc(t,:) = [lc, sc, sh, df];
    PST_nc(t,:) = dskt(YY, lc, sc, sh, df);
    QST_nc(t,:) = qskt(QQ, lc, sc, sh, df);
    PS_nc(t) = dskt(OT_nc(t), lc, sc, sh, df);
    
    qqTemp = qskt([delta1:delta1:alpha], Par_nc(t,1), Par_nc(t,2), Par_nc(t,3), Par_nc(t,4));
    SF_nc(t,:) = 1/alpha * sum(qqTemp * delta1);
    
    qqTemp = qskt([1-delta1:-delta1:1-alpha], Par_nc(t,1), Par_nc(t,2), Par_nc(t,3), Par_nc(t,4));
    LR_nc(t,:) = 1/alpha * sum(qqTemp * delta1);
    
    % Forecast
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_fc(t, :), QQ);
    Par_fc(t,:) = [lc, sc, sh, df];
    PST_fc(t,:) = dskt(YY, lc, sc, sh, df);
    QST_fc(t,:) = qskt(QQ, lc, sc, sh, df);
    PS_fc(t) = dskt(OT_fc(t), lc, sc, sh, df);
    
    qqTemp = qskt([delta1:delta1:alpha], Par_fc(t,1), Par_fc(t,2), Par_fc(t,3), Par_fc(t,4));
    SF_fc(t,:) = 1/alpha * sum(qqTemp * delta1);
    
    qqTemp = qskt([1-delta1:-delta1:1-alpha], Par_fc(t,1), Par_fc(t,2), Par_fc(t,3), Par_fc(t,4));
    LR_fc(t,:) = 1/alpha * sum(qqTemp * delta1);
    
end

save("output/matlab_output/output_" + mode)