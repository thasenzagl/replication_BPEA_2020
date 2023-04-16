function [F,Lambda,eigvals,rss,res] = pca(X,r)

    % Principal components analysis
    
    T = size(X,1);
    [V,D] = eigs(X*X',r);
    F = sqrt(T)*V; % Factors
    Lambda = (X'*F)/T; % Loadings
    eigvals = diag(D); % Eigenvalues
    
    res = X-F*Lambda'; % Residuals
    rss = res(:)'*res(:); % RSS
    
end