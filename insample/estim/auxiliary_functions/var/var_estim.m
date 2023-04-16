function [Ahat, Sigmahat, XpX, T_eff] = var_estim(Y, p)

    % VAR(p) least-squares estimator, no intercept

    T_eff = size(Y,1)-p; % Effective sample size
    X = lagmatrix(Y,1:p);
    XpX = X(p+1:end,:)'*X(p+1:end,:);
    
    Ahat = (X(p+1:end,:)\Y(p+1:end,:))';
    res = Y(p+1:end,:) - X(p+1:end,:)*Ahat';
    Sigmahat = (res'*res)/T_eff; 

end