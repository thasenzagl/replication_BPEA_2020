function conddens = ComputeCGMCondDensity(X, y, gamma, delta)
%% ComputeCGMCondDensity: compute conditional density for Gaussian model
% 
% Description: ComputeCGMCondDensity evaluates the conditional
% density/sequence of densities for the conditionally Gaussian model at the
% given values of the response variable.
% 
% Input arguments:
% - X : T-by-k matrix containing values of the conditioning variables.
% - y : Vector of points to evaluate the conditional density. Can be either
%       a row or column vector (if y is a column vector, it should have the
%       same number of rows as X.
% - gamma : Vector of length k containing coefficients for the conditional
%           mean equation.
% - delta : Vector of length k containing coefficients for the conditional
%           volatility equation.
%
% Output arguments:
% - conddens : Matrix with dimensions determined by the shape of y:
%       - If y is a column vector, conddens is a T-by-1 matrix with
%         conddens(t) = the density conditional on X(t, :) evaluated at
%         y(t)
%       - If y is a row vector (of length l), conddens is a T-by-l matrix
%         with conddens(t,j) = the density conditional on X(t, :) evaluated
%         at y(j)

condmean = X * gamma;
condsigma = exp(X * delta / 2);
conddens = normpdf(y, condmean, condsigma);

end

