function [the_mean, the_var, the_skew, the_kurt] = skew_t_moments(alpha, nu)

    % Moments of skew-t(0,1,alpha,nu) distribution:
    % mean, variance, skewness, kurtosis (not excess kurtosis)
    % See notation and formulas in Azzalini & Capitanio (JRSS 2003), p. 382
    % See also https://en.wikipedia.org/wiki/Skew_normal_distribution

    % Preliminaries
    delta = alpha./sqrt(1+alpha.^2);
    mu = delta.*(sqrt(nu/pi).*gamma(0.5*(nu-1))./gamma(0.5*nu));
    
    % Mean
    the_mean = mu;
    the_mean(nu<=1,:) = nan;
    the_mean(nu==Inf,:) = delta(nu==Inf,:)*sqrt(2/pi); % Skew-normal case

    % Variance
    the_var = nu./(nu-2) - mu.^2;
    the_var(nu<=2,:) = nan;
    the_var(nu==Inf,:) = 1 - 2*delta(nu==Inf,:).^2/pi; % Skew-normal case

    % Skewness
    the_skew = mu .* (nu.*(3-delta.^2)./(nu-3) - 3*nu./(nu-2) + 2*mu.^2) ./ (the_var.^(3/2));
    the_skew(nu<=3,:) = nan;
    the_skew(nu==Inf,:) = (4-pi)/2*(the_mean(nu==Inf,:)./sqrt(the_var(nu==Inf,:))).^3; % Skew-normal case

    % Kurtosis (not excess kurtosis)
    the_kurt = (3*nu.^2./((nu-2).*(nu-4)) - 4*mu.^2.*nu.*(3-delta.^2)./(nu-3) + 6*mu.^2.*nu./(nu-2) - 3*mu.^4) ./ (the_var.^2);
    the_kurt(nu<=4,:) = nan;
    the_kurt(nu==Inf,:) = 2*(pi-3)*(the_mean(nu==Inf,:)./sqrt(the_var(nu==Inf,:))).^4 + 3; % Skew-normal case

end