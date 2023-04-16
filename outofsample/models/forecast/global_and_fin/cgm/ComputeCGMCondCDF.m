function condcdf = ComputeCGMCondCDF(X, y, gamma, delta)
%% ComputeCGMCondCDF: compute conditional CDF for Gaussian model
%
% Description: ComputeCGMCondDensity evaluates the conditional CDF/sequence
% of CDFs for the conditionally Gaussian model at the given values of the
% response variable.
% 
% Input arguments:
% - X : T-by-k matrix containing values of the conditioning variables.
% - y : Column vector of length T containing points to evaluate the
%       conditional CDF.
% - gamma : Vector of length k containing coefficients for the conditional
%           mean equation.
% - delta : Vector of length k containing coefficients for the conditional
%           volatility equation.
%
% Output arguments:
% - condcdf : Column vector of length T with condcdf(t) = the CDF
%             conditional on X(t, :) evaluated at y(t)

condmean = X * gamma;
condsigma = exp(X * delta / 2);
condcdf = normcdf(y, condmean, condsigma);
end