function quant = skew_t_quant(p, alpha, nu, varargin)

    % Quantile of skew-t(0,1,alpha,nu) distribution

    quant = fsolve(@(y) obj(y,p,alpha,nu), tinv(p,nu), varargin{:});
    % Initial guess is based on symmetric case alpha=0

end

function [res, grad] = obj(y, p, alpha, nu)

    res = skew_t_cdf(y, alpha, nu) - p; % Residual to set to zero
    grad = skew_t_pdf(y, alpha, nu); % Gradient

end