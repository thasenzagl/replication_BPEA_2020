%______________________________________________________________________
function S = SKF_lag(Y,Z,R,T,Q,A_0,P_0)
%______________________________________________________________________
% Kalman filter for stationary systems with time-varying system matrices
% and missing data.
%
% The model is        y_t   = Z * a_t + eps_t
%                     a_t+1 = T * a_t + u_t
%
%______________________________________________________________________
% INPUT
%        Y         Data                                 (nobs x n)
% OUTPUT
%        S.Am       Predicted state vector  A_t|t-1      (nobs x m)
%        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)
%        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)
%        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)
%        S.loglik   Value of likelihood function

    % Output structure & dimensions
    [n m] = size(Z);
    nobs  = size(Y,2);

    S.Am  = ones(m,nobs)*NaN;   S.Pm  = single(ones(m,m,nobs))*NaN;
    S.AmU = ones(m,nobs+1)*NaN;   S.PmU = single(ones(m,m,nobs+1))*NaN;
    S.loglik = 0;

    %______________________________________________________________________
    Au = A_0;  % A_0|0;
    Pu = P_0;  % P_0|0

    S.AmU(:,1)    = Au;
    S.PmU(:,:,1)  = single(Pu);
    %S.PmU{1}  = Pu;



    for t = 1:nobs

        % A = A_t|t-1   & P = P_t|t-1

        A   = T*Au;
        P   = T*Pu*T' + Q;
        P   =  0.5 * (P+P');

        % handling the missing data
        [y_t,Z_t,R_t,L_t] = MissData(Y(:,t),Z,R);

        if isempty(y_t)
            Au = A;
            Pu = P;

        else
            PZ  = P*Z_t';
            iF  = pinv(Z_t*PZ + R_t);
            PZF = PZ*iF;

            V   = y_t - Z_t*A;
            Au  = A  + PZF * V;
            Pu  = P  - PZF * PZ';
            Pu   =  0.5 * (Pu+Pu');
        end

        S.Am(:,t)   = A;
        S.Pm(:,:,t) = single(P);
        %%S.Pm{t} = P;

        % Au = A_t|t   & Pu = P_t|t

        S.AmU(:,t+1)    = Au;
        S.PmU(:,:,t+1)  = Pu;
        %%S.PmU{t+1}  = Pu;
        S.loglik = S.loglik + 0.5*(log(det(iF))  - V'*iF*V);
    end % t

    if isempty(y_t)
        S.KZ = zeros(m,m);
    else
        S.KZ = PZF*Z_t;
    end
end