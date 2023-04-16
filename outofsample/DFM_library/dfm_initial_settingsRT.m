function [P, DatesV, DD, DD_series, name_vintage, StartVintage, EndVintage, blocks, target_last, include]  = dfm_initial_settingsRT(P, StartDate, EndDate)

%% Read the Dataset into Matlab using the READDATA function

load('data/real_time_data.mat')
X_Last=DD{end};

ReleasesMatVec = release_dates;
block_table = readtable(P.Path, 'Sheet', 'blocks', 'ReadVariableNames', 0);
block_cell = table2cell(block_table(2:end, 2:end));
block_cell(any(cellfun(@isempty,block_cell), 2), :) = [];
blocks = cellfun(@str2num, block_cell(:, 2:end));
include = cellfun(@str2num, block_cell(:, 1));
DatesV = datevec(ref_months);

P.all_series = DD_series{end};

iSer = find(ismember(P.all_series,P.SerNews))';
target_last = X_Last(:,iSer);

for i=1:size(DD,2)
    tM = length(DD{end}) - size(DD{i},1);
    DD{i} = [DD{i}; nan(tM, size(DD{i},2))];
end

%% Model Parameters
%
% P.r:      number of factors
% P.nM:     number of variables with monthly frequency
% P.nQ:     number of variables with quarterly frequency
% P.Series: names of the individual series

P.r=P.nR.*ones(size(blocks(1,:)));

Res_old    = [];
SubS       = [];

%% First and Last Release Dates for Vintage

start  = find(ReleasesMatVec>=datenum(StartDate),1,'first'); % first data release after StartDate
finish = find(ReleasesMatVec<=datenum(EndDate),1,'last'); % last data release before EndDate

name_vintage = ReleasesMatVec(start:finish);

%% Restrictions

P.Rconstr = [2 -1 0 0 0; ...
    3 0 -1 0 0; ...
    2 0 0 -1 0; ...
    1 0 0 0 -1];

P.Rconstr_vec = [1, 2, 3, 2, 1];

P.q = zeros(size(P.Rconstr, 1), 1);

%% Getting Ready for the Estimation

% Start and End of Vintage
StartVintage = 1;
EndVintage = length(name_vintage);

end