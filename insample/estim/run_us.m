%clear;

% Estimate models on U.S. data
% MPM 2020-01-19


%% Settings

%settings = struct;

% Specification, model, and prior
%settings.predic = 'factor';     % Predictor variables: either 'factor' or 'var';
%settings.model = 'skewt';       % Model: either 'condhet' or 'skewt';
%settings.start_year = [];       % Start year to impose on data (set to [] if full data should be used)
if strcmp(settings.predic, 'var') && strcmp(settings.model, 'condhet')
    settings.prior = 'horseshoe'; % Prior: either 'normal' or 'horseshoe' (latter only applicable for model 'condhet')
else
    settings.prior = 'normal';
end

% Data and specification
settings.dat_file = '../data/us/output/data_export.csv'; % CSV data file
settings.y_var = 'GDP';             % LHS variable Y
settings.date_var = 'date';         % Date variable
settings.factor_file = '../data/us/factors/factors.csv'; % Factors: CSV data file
settings.factor_date_var = 'Var1';  % Factors: date variable

% Stan settings
settings.stan = struct('chains', 4, ...
                       'warmup', 5e3, ...
                       'iter', 5e3, ...
                       'thin', 1, ...
                       'verbose', true, ...
                       'refresh', 100, ...
                       'inc_warmup', true, ...
                       'seed', 202002171); % Random number seed

% Output settings
settings.output_folder = fullfile('output/stan/us', settings.predic, strcat(settings.model, num2str(settings.start_year))); % Output directory


%% Load data

% Read from CSV;
dat = readtable(settings.dat_file); % GDP data
dat.dateno = datenum(dat{:,settings.date_var}, 'ddmmmyyyy'); % Convert dates to serial numbers

if strcmp(settings.predic, 'factor')
    
    % Load factor data 
    dat_factor = readtable(settings.factor_file, 'Format', '%s %f %f'); % Factor data
    dat_factor.dateno = datenum(dat_factor{:,settings.factor_date_var}, 'mm/dd/yy'); % Convert dates to serial numbers
    dat_factor = removevars(dat_factor, settings.factor_date_var);
    dat_factor_vars = dat_factor.Properties.VariableNames; % Variable names

    % Combine GDP data and factor data
    dat = innerjoin(dat, dat_factor, 'Keys', 'dateno');
    
    % Keep only selected variables
    dat = dat(:,[{settings.date_var, settings.y_var}, dat_factor_vars]);
    
end

% Restrict sample if desired
if isempty(settings.start_year)
    sample = ~isnan(dat.dateno);
else
    sample = dat.dateno>=datenum(settings.start_year,1,1);
end
dat = dat(sample,:);
dat = removevars(dat, 'dateno');


%% Run estimation

estim_model;
