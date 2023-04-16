%% Clear workspace; 
clear
close all
clc

%%
% This script runs a forecasting exercise using quantile regressions of 
% real GDP growth on the global, financial, and non-financial factors. 

% The code is nearly entirely taken from the replication 
% codes from Adrien et al. (2019).

%% Set forecast horizon (run script separately for h = 1 and h=4);
h=1;

%% Run the exercise with the global factor and financial factors
run('global_and_fin/VulnerabilityOutOfSample.m');
clearvars -except h

%% Run the exercise with the global factor only
run('global_only/VulnerabilityOutOfSample.m');