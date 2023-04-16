function condmean = skew_t_condmean(cutoff, alpha, nu, varargin)

    % E[Y|Y<cutoff] for Y ~ skew-t(0,1,alpha,nu)

    prob = skew_t_cdf(cutoff, alpha, nu); % P(Y<cutoff)
    integr = integral(@(y) y.*skew_t_pdf(y, alpha, nu), -Inf, cutoff, varargin{:}); % E[Y*1(Y<cutoff)]
    condmean = integr/prob; % E[Y|Y<cutoff]

end