function [YQ, xbw, ybw] = ComputeNonparCondQuantiles(Xfit, yfit, QQ, Xeval, xbw, ybw)
% ComputeNonparCondCDFquantiles: nonparametric conditional quantile estimation
% 
% Description: ComputeNonparCondCDFquantiles computes nonparametric
% conditional quantiles estimates, obtained by inverting an estimated
% cumulative distribution function as in Li, Lin, and Racine (2013). If not
% provided, optimal bandwidths for the response and conditioning variables
% will be computed using the least-squares cross validation described in
% the same reference. Quantiles can be estimated for either the sample used
% to fit the conditional CDF, or for some other values provided by the
% user.
%
% Input arguments:
% - Xfit : n-by-q matrix containing values of the conditioning variables
%          used to estimate the conditional CDF. If Xeval is not provided,
%          the conditional quantiles returned by the function will be
%          computed for these values of the conditioning variables.
% - yfit : Vector of length n containing values of the response variable
%          used to estimate the conditional CDF.
% - QQ : Vector of numbers between 0 and 1 (exclusive) indicating the
%        quantiles that will be estimated. Defaults to 0.05:0.05:0.95.
% - Xeval : neval-by-q matrix containing values of the conditioning 
%           variables for which the conditional quantiles will be
%           evaluated. Defaults to Xfit (i.e. in-sample estimation).
% - ybw : Positive float, bandwidth corresponding to y used to compute the
%         nonparametric CDF estimate. Defaults to the optimal values
%         obtained by minimizing the least-squares cross-validation
%         function of Li, Lin, and Racine (2013) using yfit and Xfit.
% - xbw : Vector of length n containing positive floats, bandwidths
%         corresponding to y used to compute the nonparametric CDF
%         estimate. Defaults to the optimal values obtained by minimizing
%         the least-squares cross-validation function of Li, Lin, and
%         Racine (2013) using yfit and Xfit.
%
% Output arguments:
% - YQ : neval-by-length(QQ) matrix containing the estimated conditional
%        quantiles, evaluated for the values of the conditioning variables
%        in each row of Xeval.
% - ybw : positive float, the bandwidth corresponding to y used to compute
%         the estimated quantiles.
% - xbw : Vector of length n containing positive floats, the bandwidths
%         corresponding to x used to compute the estimated quantiles.

%% Settings
if nargin < 3
    QQ = 0.05:0.05:0.95;
end
if nargin < 4
    Xeval = Xfit;
end
if nargin < 6
    % Manually compute bandwidths if both ybw, xbw not provided
    [xbw, ybw] = ComputeNonparCondCDFbw(Xfit, yfit, QQ);
end

neval = size(Xeval,1);

[n, q] = size(Xfit);
numQQ = length(QQ);
YQ = NaN(neval, numQQ);

%% To speed up iterations, precompute differences and kernel values for x
% xdiff(i, j, s) = Xeval(i, s) - Xfit(j, s)
xdiff = NaN(neval, n, q);
for i = 1:neval
    for j = 1:n
        xdiff(i, j, :) = Xeval(i, :) - Xfit(j, :);
    end
end
% xbwmat(i, j, s) = xbw(s)
xbwmat(1, 1, :) = xbw; xbwmat = repmat(xbwmat, neval, n, 1);
% Kmat(i, j) = K_h(Xeval(i, :), Xfit(j, :))
Kmat = prod(normpdf(xdiff ./ xbwmat), 3) / prod(xbw);
% sumKmat = (unscaled) kernel density estimate of marginal density of x
sumKmat = sum(Kmat, 2);

%% Compute conditional quantiles for each x_i
% Initial guess for smallest conditional quantile
% (unconditional quantile of y)
yq0 = quantile(yfit, QQ(1));

% Compute F^-1 (QQ(qind) | x_i) for each observation and quantile
for i = 1:neval
    yq = yq0;
    for qind = 1:numQQ
        qq = QQ(qind);        
        
        % Create inline function for conditional density minus quantile
        fun = @(yeval) ComputeNonparCondCDF(yeval, yfit, ybw, Kmat(i, :), sumKmat(i)) - qq;
        
        yq = fzero(fun, yq); % use last quantile as initial guess for current quantile
        YQ(i, qind) = yq;
    end
end

end

function Feval = ComputeNonparCondCDF(yeval, Y, ybw, K, sumK)
% G(j) = Phi((yeval - y_j) / ybw)
G = normcdf((yeval - Y) ./ ybw);
Feval = dot(G, K) / sumK;
end