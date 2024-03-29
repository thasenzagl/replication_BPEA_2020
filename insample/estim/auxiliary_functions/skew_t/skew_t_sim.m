function T = skew_t_sim(alphas,nus)

    % Simulate from the skew-t(0,1,alpha,nu) distribution
    
    % "Distributions Generated by Perturbation of Symmetry with Emphasis on
    % a Multivariate Skew t-Distribution"
	% Azzalini & Capitanio (2003)
    
    % "A Probabilistic Representation of the 'Skew-Normal' Distribution"
    % Henze (1986)

    n = length(alphas);
    if length(nus)==1
        nus = nus*ones(n,1);
    end

    % Simulate skew normal
    norms = randn(n,2);
    Z = (alphas.*abs(norms(:,1))+norms(:,2))./sqrt(1+alphas.^2);
    
    % Simulate skew-t
    V = chi2rnd(nus)./nus;
    V(nus==Inf) = 1; % Skew-normal case
    T = Z./sqrt(V);

end