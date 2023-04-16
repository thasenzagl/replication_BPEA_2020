%clear;

% Estimate models on data for several OECD countries
% MPM 2020-01-20


%% Settings

%settings = struct;

% Specification, model, and prior
%settings.predic = 'factor';     % Predictor variables: either 'factor' or 'var';
%settings.model = 'condhet';       % Model: either 'condhet' or 'skewt';
settings.start_year = [];       % Start year to impose on data (set to [] if full data should be used)
if strcmp(settings.predic, 'var') && strcmp(settings.model, 'condhet')
    settings.prior = 'horseshoe'; % Prior: either 'normal' or 'horseshoe' (latter only applicable for model 'condhet')
else
    settings.prior = 'normal';
end

% Data and specification
settings.countries = {'AUS', 'BEL', 'CAN', 'CHE', 'DEU', 'ESP', 'FRA', 'GBR', 'ITA', 'JPN', 'NLD', 'SWE', 'USA'}; % Country codes
settings.dat_folder = '../data/global/output/csv';  % CSV data folder
settings.y_var = 'GDP';                             % LHS variable Y
settings.date_var = 'date';                         % Date variable
settings.factor_dat_folder = '../../international_factors/output/csv'; % Factors: CSV data folder
settings.factor_date_var = 'date';                  % Factors: date variable
settings.factor_dateadd = 693960;                   % Factors: Number to add to date serial numbers

% Missing data imputation
settings.r_imp = 8;         % Number of principal components used for imputation
settings.tol_imp = 1e-6;    % Convergence tolerance for imputation

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
settings.output_root_folder = fullfile('output/stan/global', settings.predic, strcat(settings.model, num2str(settings.start_year))); % Output directory


%% Run analysis by country

for i=1:length(settings.countries) % For each country...
    
    the_country = settings.countries{i};
    disp(the_country);
    
    settings.output_folder = fullfile(settings.output_root_folder, the_country);

    % Load data
    dat = readtable(fullfile(settings.dat_folder, strcat(the_country, '.csv'))); % Read country data from CSV;
    dat.dateno = datenum(dat{:,settings.date_var}, 'ddmmmyyyy'); % Convert dates to serial numbers
    
    if strcmp(settings.predic, 'factor')
        
        % Load factor data
        dat_factor = readtable(fullfile(settings.factor_dat_folder, strcat(the_country, '.csv'))); % Factor data
        dat_factor.dateno = dat_factor{:,settings.factor_date_var} + settings.factor_dateadd; % Dates (serial numbers)
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

    % Run estimation
    estim_model;

end
