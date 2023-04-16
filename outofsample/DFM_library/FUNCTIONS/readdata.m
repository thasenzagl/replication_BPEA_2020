function [X_Last, iSer, DatesV, P, Dates, ReleasesMat, ReleasesMatVec, Period, iSer_idx] = readdata(P)
%% READDATA: A Function to Read Data from .txt Files
%
%   The READDATA function takes input arguments of type structure P
%   P:  The fields of P used as inputs by READDATA are as follows:
%       StartEst is the start date of the estimation stored in a
%           vector [year month]
%       Path is the path to the folder with the text files that need to be
%           imported
%
%   The READDATA function returns the following output:
%       T: number of monthly observations starting at StartEst
%       t: number of monthly observations starting in the first
%          month with data in the Calander worksheet
%   X_Last:         Txn matrix of transformed monthly and quarterly data
%   SeriesID:       nx1 vector of IDs of included monthly and quarterly series
%   Group:          nx1 of the group of included monthly series
%   DatesV:         Tx2 matrix of dates starting at in the month
%                       initialized by P.StartEst. The first column is the year, the second
%                       column is the month
%   nQ:             Dimensions of transformed, quarterly data matrix
%   nM:             Dimensions of transformed, monthly data matrix
%   P:              Stucture containing fields used as inputs and nx3 matrix
%                       P.blocks. If an entry is 1 the corresponding series
%                       is part of the block global, hard or soft, respectively
%   Dates:           Tx1 matrix of monthly dates saved as Matlab serial date numbers starting at in the month
%                       initialized by P.StartEst
%   Descr:          nx1 vector with descriptions of the variables
%   Info:           Contains information on the series such as Release Unit
%                   and Release Frequency
%   Included:       nx1 vector that indicates if monthly variables are included in model
%   ReleasesMat:      nxt matrix. Like ReleasesM but in Matlab date format.
%   ReleasesMatVec: Vector of release dates. The number of entries is give by
%                       the number of unique release dates. Sorted and does
%                       not contain duplicates. In Matlab date format.
%   Period:         tx6 matrix of dates in which data was released. Beginns at
%                       StartDate and ends at EndDate.
%   Trform:         nx3 matrix that shows which transformations are
%                   performed for the individual series

    %% Initialize Input Variables
    %
    StartEst = P.StartEst;
    path = P.Path;

    %% Loading Monthly Data: Worksheet block_m
    %
    %   IncludedM   Indicates if monthly variables are included in model
    %   BlocksM:    Column 1 is 1 if series is global data, 0 otherwise
    %               Column 2 is 1 if series is hard data, 0 otherwise
    %               Column 3 is 1 if series is soft data, 0 otherwise
    %   SeriesIDM:  IDs of included monthly series
    %   GroupM:     The group of included monthly series

    % Read in Data    
    block_m = readtable(path, 'Sheet', ...
                        'block_m', 'ReadVariableNames', 0);
    
    block_m_cell = table2cell(block_m(2:end, 2:end));
    block_m_cell(any(cellfun(@isempty,block_m_cell), 2), :) = [];

    % Variables
    IncludedM = cellfun(@str2num, block_m_cell(:, 1));
    BlocksM = cellfun(@str2num, block_m_cell(:, 2:end));
    BlocksM = BlocksM(IncludedM == 1,:);
    SeriesIDM = table2cell(block_m(2:end, 1));
    SeriesIDM = SeriesIDM(IncludedM == 1,:);

    clearvars block_m_cell

    
    %% Loading Monthly Data: Worksheet LegendMonthly
    %
    %   TrforM: Column 1 is 1 if series was log transformed, 0 otherwise
    %                  Column 2 is 1 if first difference of series were taken, 0 otherwise
    %                  Column 3 is 1 if series was filtered, 0 otherwise

    % Read in Data
    LegendMonthly = readtable(path, 'Sheet', ...
                              'LegendMonthly', 'ReadVariableNames', 0);

    % Variables
    TrformM = table2cell(LegendMonthly(2:end, 2));
    TrformM = TrformM(IncludedM == 1);
    TrformM = cellfun(@str2num, TrformM);

    %% Loading Monthly Data: Worksheet monthly
    %
    %   DataM:      Matrix of monthly data starting at in the month
    %                   initialized by P.StartEst
    %   DatesM:     Monthly dates saved in a vector of Matlab serial date numbers starting at in the month
    %                   initialized by P.StartEst
    %   DatesMV:    Monthly dates saved in a date vector starting at in the month
    %                   initialized by P.StartEst
    %   T:          Number of months starting at StartEst

    % Read in Data    
    monthly = readtable(path, 'Sheet', ...
                        'monthly', 'ReadVariableNames', 0);
                            
    monthly_cell = table2cell(monthly(2:end , :));
    monthly_cell(any(cellfun(@isempty,monthly_cell(:, 1)), 2), :) = [];
    
    % Variables
    DatesM = cell2mat(monthly_cell(:, 1)) + datenum('12-30-1899');  % Excel date number to Matlab date number 
    DatesMV = datevec(DatesM);
    
    DataM = monthly_cell(:, 2:end);
    empties = cellfun('isempty',DataM);
    DataM(empties) = {'NaN'};
    DataM = cellfun(@str2num, DataM);
    DataM = DataM(:, IncludedM == 1);
    T = length(DatesM);
    
    clearvars monthly_cell monthly_cell_info
    
    %% Loading Quarterly Data: Worksheet block_q
    %
    %   IncludedQ:  Indices of quarterly variables that are included in model
    %   BlocksQ:    Column 1 is 1 if series is global data, 0 otherwise
    %               Column 2 is 1 if series is hard data, 0 otherwise
    %               Column 3 is 1 if series is soft data, 0 otherwise
    %   SeriesIDQ:  IDs of included monthly series

    % Initialize variables
    
    block_q = readtable(path, 'Sheet', ...
                        'block_q', 'ReadVariableNames', 0);

    block_q_cell = table2cell(block_q(2:end, 2:end));
    block_q_cell(any(cellfun(@isempty,block_q_cell),2),:) = [];

    % Variables
    IncludedQ = cellfun(@str2num, block_q_cell(:, 1));
    BlocksQ  = block_q_cell(:, 2:end);
    BlocksQ  = cellfun(@str2num, BlocksQ(IncludedQ == 1,:));
    SeriesIDQ = table2cell(block_q(2:end, 1));
    SeriesIDQ = SeriesIDQ(IncludedQ == 1,:);
    
    clearvars block_q_cell

    
    %% Loading Quarterly Data: Worksheet LegendQuarterly

    % Read in Data
    LegendQuarterly = readtable(path, 'Sheet', ...
                                'LegendQuarterly', 'ReadVariableNames', 0);
            
    
    % Variables
    TrformQ = table2cell(LegendQuarterly(2:end, 2));
    TrformQ = TrformQ(IncludedQ == 1);
    TrformQ = cellfun(@str2num, TrformQ);

    %% Loading Quarterly Data: Worksheet quarterly
    %
    %   DataQ:  Quartlery released data that is included in the model
    %

    % Read in Data
    quarterly = readtable(path, 'Sheet', ...
                          'quarterly', 'ReadVariableNames', 0);

    % Variables (Quarterly frequency)
    quarterly_cell = table2cell(quarterly(2:end, :));
    
    DataQ = quarterly_cell(:, 2:end);
    empties = cellfun('isempty',DataQ);
    DataQ(empties) = {'NaN'};
    DataQ = cellfun(@str2num, DataQ);
    DataQ = DataQ(:, IncludedQ == 1);
    
    clear quarterly_cell quarterly_cell_info

    %% Loading Calendar (Release Dates)
    %
    %   Included:       Indices of monthly and quarterly variables that are
    %                       included in model.
    %   ReleasesE:      nxt matrix. Every column is a vector of the release dates of the
    %                       included data series from one specific month.
    %   ReleasesM:      Like ReleasesE but in Matlab format.
    %   ReleasesMVec:   Vector of release dates. Sorted and does not contain
    %                       duplicates. In Matlab date format.
    %

    % Read in Data
    calendar = readtable(path, 'Sheet', ...
                          'calendar', 'ReadVariableNames', 0);

    % Temporary variables
    
    calendar_cell = table2cell(calendar(2:end, 2:end));
    calendar_cell(any(cellfun(@isempty,calendar_cell),2),:) = [];
        
    % Final variables
    ReleasesEVec = cell2mat(calendar_cell(:));
    Included = [IncludedM;IncludedQ];
    ReleasesExcel = cell2mat(calendar_cell(Included == 1, :));
    ReleasesMat = ReleasesExcel + datenum('12-30-1899');
    ReleasesMatVec = unique(sort(ReleasesExcel(:) + datenum('12-30-1899')));

    clearvars calendar_cell calendar_cell2
    

    %% Loading Calendar (Periods)
    %
    %   Period: Matrix of dates in which data was released. Beginns at
    %               StartDate and ends at EndDate
    %   PMV:    tx2 matrix of years and months of data release dates.
    %

    % Variables
    calendar_cell3 = table2cell(calendar(1, 2:end));
    Period_temp = cell2mat(calendar_cell3(:)); % transpose row vector to get column vector
    Period = datevec(Period_temp + datenum('12-30-1899'));

    clearvars calendar_cell3

    
    %% Monthly Transformations
    %   Apply Transformations
    %
    %   DataMTrf:    Matrix of monthly, transformed data

    DataMTrf_temp = prepare_missing(DataM, TrformM);
    
    % Find dimensions of matrix DataMTrf_temp
    [tM,nM] = size(DataMTrf_temp);

    % Final data matrices. If transformed matrix is has fewer rows than
    % untransformed data, add rows of NaN.
    DataMTrf = [DataMTrf_temp; NaN(T-tM,nM)];

    clearvars DataMTrf_tempDataQTrf_temp

    
    %% Quarterly Transformations
    %   Apply Transformations
    %
    %   DataQTrf:    Matrix of quarterly, transformed data

    DataQTrf_temp = prepare_missing(DataQ, TrformQ);

    % Quarterly to Monthly for DataQTrf
    DataQTrf_temp = kron(DataQTrf_temp,[NaN;NaN;1]);
    
    % Find dimensions of matrix DataMTrf_temp
    [tQ,nQ] = size(DataQTrf_temp);
    
    % Final data matrices. If transformed matrix is has fewer rows than
    % untransformed data, add rows of NaN.
    DataQTrf = [DataQTrf_temp; NaN(T-tQ,nQ)];
    
    clearvars DataQTrf_temp
    
    %% Final Modifications
    %
    %   indexEst:   Index of date vector that corresponds to StartEst

    indexEst = find(DatesMV(:,1) == StartEst(1) & DatesMV(:,2) == StartEst(2));

    % Combining Data Matrices for monthly and quarterly
    Data_temp = [DataMTrf DataQTrf];
    Data_temp = real(Data_temp);
    SeriesID = [SeriesIDM; SeriesIDQ];
    P.blocks = [BlocksM; BlocksQ];
    
    % Series start at time StartEst
    X_Last = Data_temp(indexEst:end,:);
    DatesV = DatesMV(indexEst:end,1:2);
    Dates = DatesM(indexEst:end,:);

    clear Data_temp

    % Quarterly to Monthly for DataQ
    DataQ = kron(DataQ,[NaN;NaN;1]);

    % Find dimensions of matrix DataQ
    [tQ2,nQ2] = size(DataQ);
    DataQ = [DataQ; NaN(T-tQ2,nQ2)];

    % Find dimensions of matrix DataM
    [tM2,nM2] = size(DataM);
    DataM = [DataM; NaN(T-tM2,nM2)];

    %% Remove target
    
    P.Series=SeriesID;

    % The index of the variable we are nowcasting (GDP)
    iSer_idx = ~ismember(P.Series,P.SerNews);
    iSer = find(ismember(P.Series,P.SerNews))';
        
    % nM and nQ
    if iSer>nM
        nQ = nQ - size(iSer,1);
    else
        nM = nM - size(iSer,1);
    end    
    
    P.nQ = nQ;
    P.nM = nM;
    
    % Blocks and Series
    P.blocks = P.blocks(iSer_idx,:);
    P.Series = P.Series(iSer_idx,:);
    
end