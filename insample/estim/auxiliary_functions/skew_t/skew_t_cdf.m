function cdf = skew_t_cdf(y, alpha, nu)

    % Distribution function of skew-t(0,1,alpha,nu) distribution, evaluated at y
    % See notation and formulas in Azzalini & Capitanio (JRSS 2003), p. 381
    
    % Preliminaries
    delta = alpha/sqrt(1+alpha^2);
    C = [1 -delta; -delta 1];
    X = [zeros(size(y(:))), y(:)];
    
    % Evaluate CDF using multivariate Student-t CDF
    if nu==Inf
        cdf = 2*mvncdf(X, [], C);
    else
        cdf = 2*mvtcdf(X, C, nu);
    end

end