function [vint, DatesV, iSer, iSer_idx, P] = readdata_fred(P, mode)

disp("Loading real time vintage for "+ P.name_vintage)

addpath("FredFetch");

%% Read excel file

% load the series IDs and transformations from the excel file
series_sheet = readtable(P.Path, 'Sheet', 'series', 'ReadVariableNames', 1);
haver_data_sheet = readtable(P.Path, 'Sheet', 'haver_data', 'ReadVariableNames', 1);
haver_dates_sheet = readtable(P.Path, 'Sheet', 'haver_dates', 'ReadVariableNames', 1);

dates = table2array(haver_data_sheet(:,1)) + datenum('12-30-1899');
idx(1) = find(dates <= datenum(P.start_date), 1,'first');
idx(2) = find(dates >= datenum(P.start_date), 1,'last');
dates = dates(idx(1):idx(2));
calendar = table2array(haver_dates_sheet) + datenum('12-30-1899');

SeriesID = table2cell(series_sheet(:,1));
transform = table2array(series_sheet(:,2));

if strcmp(mode,'nonfin_only')
    P.blocks = ones(sum(table2array(series_sheet(:,5))),1);
else    
    P.blocks = table2array(series_sheet(:,3:4));
end

%% Load vintages from Fred
vint = nan(size(dates,1), size(SeriesID,1));

for i = 1:size(SeriesID,1)
    
    
    if ismember(SeriesID{i}, haver_dates_sheet.Properties.VariableNames)
            series_idx = find(strcmp(SeriesID{i}, haver_dates_sheet.Properties.VariableNames));
            value = table2array(haver_data_sheet(:,series_idx));
               
            series_transf = prepare_missing(value, transform(i));
            ref_idx = find(calendar(:,series_idx)<=datenum(P.name_vintage), 1,'last');
            
            idx(1) = 1;
            idx(2) = find(dates==calendar(ref_idx,1));
            vint(idx(1):idx(2),i) = series_transf(idx(1):idx(2));
               
    else
        s = fred.vint(SeriesID{i}, P.name_vintage);
        
        if strcmp(s.frequency,'M')
            
            idx1(1) = find(s.date>=datenum(P.start_date),1,'first');
            idx1(2) = find(s.date<=datenum(P.end_date),1,'last');
            idx2(1) = find(dates==s.date(idx1(1)));
            idx2(2) = find(dates==s.date(idx1(2)));
            
            series_transf = prepare_missing(s.value, transform(i));
            vint(idx2(1):idx2(2),i) = series_transf(idx1(1):idx1(2));
            
        elseif strcmp(s.frequency,'Q')
            date = kron(s.date,[NaN;NaN;1]);
            
            idx1(1) = find(date>=datenum(P.start_date),1,'first');
            idx1(2) = find(date<=datenum(P.end_date),1,'last');
            idx2(1) = find(dates==date(idx1(1)));
            idx2(2) = find(dates==date(idx1(2)));
            
            series_transf = prepare_missing(s.value, transform(i));
            series_transf = kron(series_transf,[NaN;NaN;1]);
            
            vint(idx2(1)+2:idx2(2)+2,i) = series_transf(idx1(1):idx1(2));
        end
    end
end

%% Dates
DatesV = datevec(dates);
DatesV = DatesV(:,1:2);

%% Select data for nonfinancial factor
if strcmp(mode,'nonfin_only')
    idx=logical(table2array(series_sheet(:,5)));
    SeriesID=SeriesID(idx);
    vint = vint(:,idx);
end    

%% Remove target

P.Series=SeriesID;

% The index of the variable we are nowcasting (GDP)
iSer_idx = ~ismember(P.Series,P.SerNews);
iSer = find(ismember(P.Series,P.SerNews))';

P.nQ = 0;
P.nM = sum(iSer_idx);

% Blocks and Series
P.blocks = P.blocks(iSer_idx,:);
P.Series = P.Series(iSer_idx,:);

%% Estimation parameters and restrictions

P.r=P.nR.*ones(size(P.blocks(1,:)));

P.Rconstr = [2 -1 0 0 0; ...
    3 0 -1 0 0; ...
    2 0 0 -1 0; ...
    1 0 0 0 -1];

P.Rconstr_vec = [1, 2, 3, 2, 1];

P.q = zeros(size(P.Rconstr, 1), 1);

% A vector of logical 1 (true) for monthly variables and logical 0 for
P.i_idio = logical([ones(P.nM,1);zeros(P.nQ,1)]);
end