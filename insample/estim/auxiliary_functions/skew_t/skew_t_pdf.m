function pdf = skew_t_pdf(y, alpha, nu)

    % Density of skew-t(0,1,alpha,nu) distribution, evaluated at y
    % See notation and formulas in Azzalini & Capitanio (JRSS 2003), p. 380
    
    if nu==Inf
        pdf = 2 * normpdf(y) .* normcdf(alpha.*y); % Skew-normal case
    else
        pdf = 2 * tpdf(y,nu) .* tcdf(alpha.*y.*sqrt((nu+1)./(y.^2+nu)), nu+1);
    end

end