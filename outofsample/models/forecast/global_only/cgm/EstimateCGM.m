function [gammahat, deltahat, loglike, hessian] = EstimateCGM(X, y, gammainit, deltainit)
%% EstimateCGM: estimate Gaussian model via maximum likelihood
% 
% Description: EstimateCGM computes the maximum likelihood estimates of the
% parameters for the conditionally Gaussian model with conditioning
% variables in both the mean and volatility equations, as described in
% Section 4.2.
%
% Input arguments:
% - X : T-by-k matrix containing values of the conditioning variables.
% - y : Vector of length T containing values of the response variable.
% - gammainit : Optional vector of length k, initial value for the
%               parameters of the mean equation, to use for numerical
%               maximization of the likelihood function. Defaults to the
%               OLS estimate of gamma.
% - deltainit : Optional vector of length k, initial value for the
%               parameters of the volatility equation, to use for numerical
%               maximization of the likelihood function. Defaults to a
%               vector of zeros.
% 
% Output arguments:
% - gammahat : Vector of length k, maximum likelihood estimate for the
%              parameters of the mean equation.
% - deltahat : Vector of length k, maximum likelihood estimate for the
%              parameters of the volatility equation.
% - loglike : Float, value of the negative log likelihood function
%             evaluated at the maximum likelihood estimates of the
%             parameters.
% - hessian : 2k-by-2k matrix, Hessian of the negative log likelihood
%             function evaluated at the maximum likelihood estimates of the
%             parameters.

k = size(X, 2);

% default initial parameter estimates
if nargin < 4
    % OLS estimate for gamma 
    gammainit = (X' * X) \ X' * y;
    % assuming intercept term is included: set intercept to log variance of
    % OLS residuals, other coefficients to zero
    deltainit = zeros(k, 1); deltainit(1) = log(var(y - X * gammainit));
end

objFun = @(par) -1 * sum(log(ComputeCGMCondDensity(X, y, par(1:k), par((k+1):end))));
par0 = [gammainit; deltainit];

[par, loglike, ~, ~, ~, hessian] = fminunc(objFun, par0);
gammahat = par(1:k);
deltahat = par((k+1):end);
loglike = -1 * loglike;
% MLE standard errors can be computed as sqrt(diag(inv(hessian)))
end
