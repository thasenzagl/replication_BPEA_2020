function [xbw, ybw, CV] = ComputeNonparCondCDFbw(X, y, xbw0, ybw0, nmulti)
%% ComputeNonparCondCDFbw: compute nonparametric conditional CDF bandwidths
% 
% Description: ComputeNonparCondCDFbw computes optimal bandwidths for
% nonparametric conditional cumulative distribution function estimation, by
% minimizing the least squares cross-validation function of Li, Lin, and
% Racine (2013). The discrete sum form of the cross-validation function is
% used, and multistarting is used for the optimization to avoid getting
% stuck at local minima.
% 
% Input arguments:
% - X : n-by-q matrix containing values of the conditioning variables.
% - y : Vector of length n containing values of the response variable.
% - xbw0 : Optional, vector of length n containing positive floats to use
%          as an initial value in the optimization for the bandwidths
%          corresponding to x. Defaults to a version of the rule-of-thumb
%          reference value.
% - ybw0 : Optional, positive float to use as an initial value in the
%          optimization for the bandwidth corresponding to y. Defaults to
%          a version of the rule-of-thumb reference value.
% - nmulti : Positive integer denoting the number of different initial
%            values to use in the multistart search for the optimal
%            bandwidths. One deterministic value initial value will be used
%            (either the user-supplied value or the rule-of-thumb reference
%            value) and the remaining nmulti-1 values will be randomly
%            chosen by MultiStart.
%
% Output arguments:
% - ybw : Positive float, the optimal bandwidth corresponding to y.
% - xbw : Vector of length n containing positive floats, the optimal
%         bandwidths corresponding to x.
% - CV : Positive float, the value of the least-squares cross-validation
%        function evaluated at the optimal bandwidths (ybw, xbw).

%% Settings
if nargin < 5
    nmulti = 5; 
end

[n, q] = size(X);

%% To speed up iterations, precompute differences for y and x
% ydiff(i, j) = y_i - y_j
ydiff = repmat(y, 1, n) - repmat(y', n, 1);
% xdiff(i, j, s) = x_i,s - x_j,s
xdiff = NaN(n, n, q);
for i = 1:n
    for j = 1:n
        xdiff(i, j, :) = X(i, :) - X(j, :);
    end
end
% Imat(i, j) = I(y_i <= y_j)
Imat = (ydiff <= 0);

%% MultiStart search for optimal bandwidth values
% Initial bandwidth values (either user-provided or rule-of-thumb default)
if nargin < 4
    yxbw0 = 1.06 * std([y, X]) * n^(-0.2 * q);
else
    yxbw0 = [ybw0, xbw0];
end

CVfun = @(yxbw) ComputeCV(yxbw(1), yxbw(2:end), ydiff, xdiff, Imat);

ms = MultiStart('Display', 'iter');
problem = createOptimProblem('fminunc', 'objective', CVfun, 'x0', yxbw0);
[yxbw, CV] = run(ms, problem, nmulti);

% NOTE: ComputeCV deals with negative bandwidths by taking absolute values,
% so be sure to do the same for the optimal bandwidths returned by
% MultiStart
ybw = abs(yxbw(1));
xbw = abs(yxbw(2:end));

% Multiply CV by n * (n-1) to conform with the original expression in Li,
% Lin and Racine (2013)
CV = CV / (n * (n-1));
end

function CV = ComputeCV(ybw, xbw, ydiff, xdiff, Imat)
%% ComputeCV: compute least-squares cross-validation function for given bws

% In case negative bandwidths are provided as arguments, take absolute
% values (and later, make sure to convert the optimal bandwidths returned
% by MultiStart to absolute values)
ybw = abs(ybw);
xbw = abs(xbw);

n = size(ydiff, 1);
logicalDiag = logical(eye(n));

% xbwmat(i, j, s) = xbw(s)
xbwmat(1, 1, :) = xbw; xbwmat = repmat(xbwmat, n, n, 1);

% Evaluate G((y_i - y_j) / h_0) and K_h(x_i, x_j) for all pairs (i, j).
% Diagonal entries of Kmat will be set = 0 for convenient computation of
% leave-one-out conditional density estimates.
Gmat = normcdf(ydiff ./ ybw);
Kmat = prod(normpdf(xdiff ./ xbwmat), 3) / prod(xbw); Kmat(logicalDiag) = 0;

% Fminusi(i, j) = Fhat_{-i} (y_j | x_i) is the leave-one-out
% conditional CDF estimator evaluated at y_j | x_i
Fminusi = (Kmat * Gmat') ./ sum(Kmat)'; %Fminusi(logicalDiag) = 1;
% since Imat(i, i) == 1, the second part ^ is a convenient way to
% exclude the (i, i) terms when computing the double sum over (i, j) for CV
% To replicate results from R, comment out the second part

% CV is not normalized by n(n-1), in order to speed up computation
CV = sum(sum((Fminusi - Imat).^2));
end