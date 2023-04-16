clear all;

% Replication files: When is Growth at Risk?
% Step 3: Create figures and tables


%% Settings

run_spec = 1; % Specification to estimate, number between 1 and 8 (see below)


%% Run estimation

start_year = [];
zero_factor = [];

switch run_spec
    
    % Figures 10-11, S.3, S.5-S.7
    case 1 % U.S. data, predictors: factors, model: skew-t
        spec = 'us';
        predic = 'factor';
        model = 'skewt';
        
    % Figure S.4
    case 2 % U.S. data, predictors: factors, model: skew-t, start year: 1980
        spec = 'us';
        predic = 'factor';
        model = 'skewt';
        start_year = 1980;
        
    % Figures S.10-S.12
	case 3 % Cross-country data, predictors: factors, model: skew-t
        spec = 'global';
        predic = 'factor';
        model = 'skewt';
        
    % Figures 12, S.13
    case 4 % U.S. data, predictors: individual variables, model: conditional heteroskedasticity
        spec = 'us';
        predic = 'var';
        model = 'condhet';
    
    % Prepare results later used to create Figures 13, S.14 and Tables 2-3
    case 5 % Cross-country data, predictors: individual variables, model: conditional heteroskedasticity
        spec = 'global';
        predic = 'var';
        model = 'condhet';
    
	% Prepare results later used to create Tables 4, S.5
    case 6 % Cross-country data, predictors: individual variables, model: skew-t
        spec = 'global';
        predic = 'var';
        model = 'skewt';
    
    % Figure S.8
    case 7 % U.S. data, predictors: factors, model: skew-t, zeroing out global factor
        spec = 'us';
        predic = 'factor';
        model = 'skewt';
        zero_factor = 1;
        
    % Figure S.9
    case 8 % U.S. data, predictors: factors, model: skew-t, zeroing out financial factor
        spec = 'us';
        predic = 'factor';
        model = 'skewt';
        zero_factor = 2;
    
end

run('estim/results_report');
