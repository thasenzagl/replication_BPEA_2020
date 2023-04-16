function [P, target_last, DatesV, DD, name_vintage, StartVintage, EndVintage, Period, ReleasesMat, iSer, iSer_idx]  = dfm_initial_settings(P, StartDate, EndDate)

%% Read the Dataset into Matlab using the READDATA function
[X_Last, iSer, DatesV, P, Dates, ReleasesMat, ReleasesMatVec, Period, iSer_idx] = readdata(P);
target_last = X_Last(:,iSer);
%% Model Parameters
%
% P.r:      number of factors
% P.nM:     number of variables with monthly frequency
% P.nQ:     number of variables with quarterly frequency
% P.Series: names of the individual series

P.r=P.nR.*ones(size(P.blocks(1,:)));

Res_old    = [];
SubS       = [];

%% First and Last Release Dates for Vintage
%
% name_vintage: A vector of all data releases between StartDate and EndDate

start  = find(ReleasesMatVec>=datenum(StartDate),1,'first')-1; % first data release after StartDate
finish = find(ReleasesMatVec<=datenum(EndDate),1,'last'); % last data release before EndDate

name_vintage = ReleasesMatVec(start:finish);


%% Create DD array and DDname matrix
%
%  DD:     Array of Matrices, one for each date at which a datapoint was
%               released (time t) in the period of a given vintage (period between start
%               and finish). Each matrix contains the monthly data of the monthly
%               and quarterly series before and including the release at time t.
%               This way, every matrix in DD is equivalent to the previous matrix
%               except that it includes the new datapoint released at date t.
%   DDname: The name of each matrix in DD is equal to the release date it
%               corresponds to.
%

% Preallocate Arrays
V = NaN(size(X_Last));
DD = cell(1, size(name_vintage,1));
DDname = cell(1, size(name_vintage,1));

% Loop over release dates and months and create DD and DDname arrays
for t=1:size(name_vintage,1);
    for i=1:size(ReleasesMat,1);
        c=find(ReleasesMat(i,:)<=name_vintage(t));
        cm = max(c);
        index = find(DatesV(:,1) == Period(cm,1) & DatesV(:,2) == Period(cm,2));
        V(1:index,i) = X_Last(1:index,i);
    end
    DD{t}= V;
    DDname{t}=name_vintage(t);
end
X_old = DD{1};
DD=DD(2:end);
name_vintage=name_vintage(2:end);
DDname=DDname(2:end);

%% Restrictions
%
% We link the unobserved monthly growth rate of the quarterly variable
% (GPD) y_t with the observed quarterly (GDP) data.
%
% y_Qt = Y_Qt-Y_Q(t-3) = [Y_Qt+Y_Q(t-1)+Y_Q(t-2)]-[Y_Q(t-3)+Y_Q(t-4)+Y_Q(t-5)]
%      = y_t + 2y_(t?1) + 3y_(t?2) + 2y_(t?3) + y_(t?4)
%
% For details, see Babura et al, 2011 (Oxford handbook chapter), Appendix B.

% -> Mariano Murasawa 2013

P.Rconstr = [2 -1 0 0 0; ...
    3 0 -1 0 0; ...
    2 0 0 -1 0; ...
    1 0 0 0 -1];

P.Rconstr_vec = [1, 2, 3, 2, 1];

P.q = zeros(size(P.Rconstr, 1), 1);

% A vector of logical 1 (true) for monthly variables and logical 0 for
P.i_idio = logical([ones(P.nM,1);zeros(P.nQ,1)]);

%% Getting Ready for the Estimation

% Start and End of Vintage
StartVintage = 1;
EndVintage = length(name_vintage);

end