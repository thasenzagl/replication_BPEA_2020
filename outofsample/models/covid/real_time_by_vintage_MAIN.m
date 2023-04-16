clear all;
close all;
clc;

%%
% This script runs a real time forecasting exercise on real GDP
% growth using quantile regressions. First we download a real-time
% vintagefrom Fred. Then we preform the following steps:

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

% Download the vitage from this date
name_vintage = {'2020-2-3', '2020-3-2', '2020-4-1',};

% GDP Fred code
P.SerNews = 'GDPC1';

% First and last date in the vintage
P.start_date = '1980-01-01';
P.end_date = '2021-12-01';

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
P.Path = "data/real_time.xlsx";

%% Quantile regression settings

deltaYY = 0.1;
YY = (-20):deltaYY:20;
QQ = 0.05:0.05:0.95;

% Find the location of the median in QQ
[~,jq50] = min(abs(QQ-.50));

%% Load the data from Fred and build the vintage
DD = cell(size(name_vintage,2),1);
YQ_tm1 = nan(size(name_vintage,2), length(QQ));
YQ_t0 = nan(size(name_vintage,2), length(QQ));
YQ_t1 = nan(size(name_vintage,2), length(QQ));
YQ_t2 = nan(size(name_vintage,2), length(QQ));
YQ_t3 = nan(size(name_vintage,2), length(QQ));
YQ_t4 = nan(size(name_vintage,2), length(QQ));

Par_tm1 = nan(size(name_vintage,2), 4);
Par_t0 = nan(size(name_vintage,2), 4);
Par_t1 = nan(size(name_vintage,2), 4);
Par_t2 = nan(size(name_vintage,2), 4);
Par_t3 = nan(size(name_vintage,2), 4);
Par_t4 = nan(size(name_vintage,2), 4);

PST_tm1 = nan(size(name_vintage,2), length(YY));
PST_t0 = nan(size(name_vintage,2), length(YY));
PST_t1 = nan(size(name_vintage,2), length(YY));
PST_t2 = nan(size(name_vintage,2), length(YY));
PST_t3 = nan(size(name_vintage,2), length(YY));
PST_t4 = nan(size(name_vintage,2), length(YY));

for t=1:size(name_vintage,2)
    P.name_vintage = name_vintage(t);
    [DD{t}, DatesV, iSer, iSer_idx, P] = readdata_fred(P, mode);
    
    if t==1
        factor_global = nan(floor(length(DatesV)/3), size(name_vintage,2));
        factor_fin = nan(floor(length(DatesV)/3), size(name_vintage,2));
    end
    
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
    X=DD{t}(:,iSer_idx);
    
    if strcmp(mode,'global_only') || strcmp(mode,'global_and_fin')
        [fg, ff, y, datesQ, jt, first_not_na] = estim_factor(X, P, DatesV, iQ, target, annualize);
        factor_global(:,t) = fg;
        factor_fin(:,t) = ff;
    elseif strcmp(mode,'nonfin_only')
        [fg, ~, y, datesQ, jt, first_not_na] = estim_factor(X, P, DatesV, iQ, target, annualize);
        factor_global(:,t) = fg;
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
            YQ_tm1(t,jq) = Z(jt-1, :) * b;
        else
            b = rq(Z(first_not_na:(jt - 1), :), y(first_not_na:(jt - 1)), QQ(jq));
        end
        
        % 0, 1, and 4 quarters ahead
        YQ_t0(t, jq) = Z(jt, :) * b;
        YQ_t1(t, jq) = Z(jt+1, :) * b;
        YQ_t2(t, jq) = Z(jt+1, :) * b;
        YQ_t3(t, jq) = Z(jt+1, :) * b;        
        YQ_t4(t, jq) = Z(jt+4, :) * b;
    end
    
    %Fit student t for tm1
    if isnan(y(jt-1))
        [lc, sc, sh, df] = QuantilesInterpolation(YQ_tm1(t, :), QQ);
        Par_tm1(t,:) = [lc, sc, sh, df];
        PST_tm1(t,:) = dskt(YY, lc, sc, sh, df);
    end
    
    %Fit student t for t0
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_t0(t, :), QQ);
    Par_t0(t,:) = [lc, sc, sh, df];
    PST_t0(t,:) = dskt(YY, lc, sc, sh, df);
    
    %Fit student t for t1
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_t1(t, :), QQ);
    Par_t1(t,:) = [lc, sc, sh, df];
    PST_t1(t,:) = dskt(YY, lc, sc, sh, df);
    
    %Fit student t for t2
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_t2(t, :), QQ);
    Par_t2(t,:) = [lc, sc, sh, df];
    PST_t2(t,:) = dskt(YY, lc, sc, sh, df);
    
    %Fit student t for t3
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_t3(t, :), QQ);
    Par_t3(t,:) = [lc, sc, sh, df];
    PST_t3(t,:) = dskt(YY, lc, sc, sh, df);
    
    %Fit student t for t4
    [lc, sc, sh, df] = QuantilesInterpolation(YQ_t4(t, :), QQ);
    Par_t4(t,:) = [lc, sc, sh, df];
    PST_t4(t,:) = dskt(YY, lc, sc, sh, df);
    
end

save("output/matlab_output/output_" + mode);
