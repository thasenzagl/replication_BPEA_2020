function YQ = ComputeCGMCondQuantiles(X, gamma, delta, QQ)
%% ComputeCGMQuantiles: compute conditional quantiles for Gaussian model
%
% Description: ComputeCGMQuantiles computes conditional quantiles for the
% conditionally Gaussian model.
% 
% Input arguments:
% - X: T-by-k matrix containing values of the conditioning variables.
% - gamma : Vector of length k containing coefficients for the conditional
%           mean equation.
% - delta : Vector of length k containing coefficients for the conditional
%           volatility equation.
% - QQ : Vector of numbers between 0 and 1 (exclusive) containing the
%        quantiles that should be computed. Defaults to 0.05:0.05:0.95
% 
% Output arguments:
% - YQ : T-by-length(QQ) matrix such that YQ(t,jq) contains the 
%        conditional quantile corresponding to QQ(jq), evaluated for the
%        values of the conditioning variables contained in X(t, :).

% Reshape vector of quantiles if given, otherwise use default quantiles.
if nargin < 4
    QQ = 0.05:0.05:0.95;
else
    QQ = reshape(QQ, 1, length(QQ));
end

% Compute conditional quantiles
cond_mean = X * gamma;
cond_sigma = exp(X * delta / 2);
for jq = 1:length(QQ)
    YQ(:,jq) = norminv(QQ(jq), cond_mean, cond_sigma);
end
end