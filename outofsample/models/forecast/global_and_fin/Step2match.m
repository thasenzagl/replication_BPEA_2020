function Res = Step2match(YQ,YQunc,QQ)
%% Step2match: Match skewed t-distributions to quantiles
%
% Description: Step2match matches skewed t-distributions to the conditonal
% quantiles provided in YQ and to the unconditional quantiles provided in
% YQunc. In addition, for each fitted distribution the PDF, CDF, moments,
% left and right entropy, and expected shortfall and longrise are computed.
% 
% The parameters of the skewed t-distribution (location, scale, shape, and
% degrees of freedom) are selected to match the following four quantiles
% provided in the input: [0.05  0.25  0.75  0.95]
%
% The grid used to compute the PDF and CDF of the fitted skewed
% t-distributions is YY = -20:0.1:20
% 
% Input arguments:
% - YQ : Matrix (with number of columns equal to length(QQ)) containing the
%        conditional quantiles to which the skewed t-distributions will be
%        fitted.
% - YQunc : Vector (with length equal to length(QQ)) containing the
%           unconditional quantiles to which the skewed t-distribution will
%           be fitted.
% - QQ : Vector of numbers between 0 and 1 (exclusive) indicating the
%        quantiles to which the columns of YQ and YQunc correspond. QQ
%        should include the quantiles listed above that will be used to
%        match the skewed t-distributions.
%
% Output arguments:
% - Res : struct containing the following fields (for convenience let T
%         denote the number of rows of YQ).
%   - PST : T-by-length(YY) matrix where PST(t,:) contains the values of
%           the probability density function of the skewed t-distribution
%           fit to the conditional quantiles YQ(t,:), evaluated at each
%           point in YY.
%   - QST : T-by-length(QQ) matrix where QST(t,:) contains the values of
%           the cumulative distribution function of the skewed
%           t-distribution fit to the conditional quantiles YQ(t,:),
%           evaluated at each quantile in QQ.
%   - STpar : T-by-4 matrix where STpar(t,:) contains the location, scale,
%             shape, and degrees of freedom parameters of the skewed-t
%             distribution fit to the conditional quantiles YQ(t,:).
%   - MeanST : Vector with length T containing the mean of the skewed
%              t-distribution fit to the conditional quantiles YQ(t,:).
%   - VarianceST : Vector with length T containing the variance of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - SkewnessST : Vector with length T containing the skewness of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - KurtosisST : Vector with length T containing the kurtosis of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - LeftEntropyST : Vector with length T containing the left entropy of
%                     the skewed t-distribution fit to the conditional
%                     quantiles YQ(t,:).
%   - RightEntropyST : Vector with length T containing the right entropy of
%                      the skewed t-distribution fit to the conditional
%                      quantiles YQ(t,:).
%   - EsfST : Vector with length T containing the expected shortfall of the
%             skewed t-distribution fit to the conditional quantiles
%             YQ(t,:).
%   - EljST : Vector with length T containing the expected longrise of the
%             skewed t-distribution fit to the conditional quantiles
%             YQ(t,:).
%   - YY : Vector containing the grid of values for which the skewed
%          t-distribution PDF and CDF are evaluated.
%   - PSTunc : Vector with length = length(YY) containing the values of
%              the probability density function of the skewed
%              t-distribution fit to the unconditional quantiles YQunc,
%              evaluated at each point in YY.
%   - QSTunc : Vector with length = length(YY) containing the values of
%              the cumulative distribution function of the skewed
%              t-distribution fit to the unconditional quantiles YQunc,
%              evaluated at each quantile in QQ.
%   - STparunc : Vector with length 4 containing the location, scale,
%                shape, and degrees of freedom parameters of the skewed-t
%                distribution fit to the unconditional quantiles YQunc.

%% Settings
% Quantiles used to match the skewed t-distribution
Qmatch = [0.05  0.25  0.75  0.95];

% Find the locations in QQ of the quantiles to be used for matching
for j=1:length(Qmatch)
    [~,J] = min(abs(QQ-Qmatch(j)));
    Jmatch(j) = J;
end
% Find the location of the median in QQ
[~,jq50] = min(abs(QQ-.50));

% Target probability for expected shortfall/longrise
alpha = 0.05;
% Interval width for approximation to expected shortfall/longrise integral
delta1 = 0.01;

% Grid of points to evaluate skewed t-density
delta = 0.1;
YY = [-20:delta:20];

%% Fit skewed t-distribution to unconditional quantiles
qqTarg = YQunc;
[lcn,scn,shn,dfn] = QuantilesInterpolation(qqTarg,QQ);

PSTunc   = dskt(YY,lcn,scn,shn,dfn);
QSTunc   = qskt(QQ,lcn,scn,shn,dfn);
STparunc = [lcn scn shn dfn];

%% Fit skewed t-distribution to conditional quantiles
% Initialize matrices to store output
PST            = NaN(size(YQ,1),length(YY));
QST            = NaN(size(YQ));
STpar          = NaN(size(YQ,1),4);
MeanST         = NaN(size(YQ,1),1);
VarianceST     = NaN(size(YQ,1),1);
SkewnessST     = NaN(size(YQ,1),1);
KurtosisST     = NaN(size(YQ,1),1);
LeftEntropyST  = NaN(size(YQ,1),1);
RightEntropyST = NaN(size(YQ,1),1);
EsfST          = NaN(size(YQ,1),1);
EljST          = NaN(size(YQ,1),1);

% Fit skewed t-distribution for each observation
for jt = 1:size(YQ,1)
    % Target quantiles, for observation jt
    qqTarg = YQ(jt,:);
    
    % Display progress message
    if mod(jt,12)==0
        disp(['Now matching quantile function for observation ',...
              num2str(jt),' (out of a total of ',num2str(size(YQ,1)),')'])
    end
    
    if ~any(isnan(qqTarg))
        
        %%% NOTE: here I no longer use previously estimated parameters as
        %%%       initial conditions
        %if jt<=5 || mod(jt,4)==0
            [lc,sc,sh,df]  = QuantilesInterpolation(qqTarg,QQ);
        %else
        %    [lc,sc,sh,df]  = QuantilesInterpolation(qqTarg,QQ,lc0,sc0,sh0,df0);
        %end
        
        %lc0=lc; sc0=sc; sh0=sh; df0=df;
        %%% END NOTE SECTION
        
        % Compute PDF and CDF for fitted distribution
        PST(jt,:)   = dskt(YY,lc,sc,sh,df);
        QST(jt,:)   = qskt(QQ,lc,sc,sh,df);
        STpar(jt,:) = [lc sc sh df];
        
        % Compute moments for fitted distribution
        MeanST(jt,:) = sum(YY .* PST(jt,:) * delta);
        VarianceST(jt,:) = sum(((YY - MeanST(jt)).^2) .* PST(jt,:) * delta);
        SkewnessST(jt,:) = sum(((YY - MeanST(jt)).^3) .* PST(jt,:) * delta) / (VarianceST(jt,:)^(3/2));
        KurtosisST(jt,:) = sum(((YY - MeanST(jt)).^4) .* PST(jt,:) * delta) / (VarianceST(jt,:)^2);
        
        % Compute left/right entropy and expected shortfall/longrise for
        % fitted distribution (using discrete sum approximation to the
        % integrals in Section 3.2)
        Temp =  PST(jt,:) .* (YY < QST(jt,jq50));
        LeftEntropyST(jt,:) = -sum((log(PSTunc) - log(PST(jt,:))) .* Temp * delta);
        
        Temp =  PST(jt,:) .* (YY > QST(jt,jq50));
        RightEntropyST(jt,:) = -sum((log(PSTunc) - log(PST(jt,:))) .* Temp * delta);
        
        qqTemp = qskt([delta1:delta1:alpha], STpar(jt,1), STpar(jt,2), STpar(jt,3), STpar(jt,4));
        EsfST(jt,:) = 1/alpha * sum(qqTemp * delta1);
        
        qqTemp = qskt([1-delta1:-delta1:1-alpha], STpar(jt,1), STpar(jt,2), STpar(jt,3), STpar(jt,4));
        EljST(jt,:) = 1/alpha * sum(qqTemp * delta1);
    end
end

%% Store all results in struct Res
Res.PST            = PST;
Res.QST            = QST;
Res.STpar          = STpar;
Res.MeanST         = MeanST;
Res.VarianceST     = VarianceST;
Res.SkewnessST     = SkewnessST;
Res.KurtosisST     = KurtosisST;
Res.LeftEntropyST  = LeftEntropyST;
Res.RightEntropyST = RightEntropyST;
Res.EsfST          = EsfST;
Res.EljST          = EljST;
Res.YY             = YY;

Res.PSTunc   = PSTunc;
Res.QSTunc   = QSTunc;
Res.STparunc = STparunc;
