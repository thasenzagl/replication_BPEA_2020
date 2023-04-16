function S = SKF_FP(Y,Z,R,T,Q,A_0,P_0)
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

    S.Am  = nan(m,nobs);   S.Pm  = nan(m,m,nobs);
    S.AmU = nan(m,nobs+1);   S.PmU = nan(m,m,nobs+1);
    S.loglik = 0;

    Plarge = nan(m,nobs);
    Tlarge = nan(m,nobs);
    
    %______________________________________________________________________
    Au = A_0;  % A_0|0;
    Pu = P_0;  % P_0|0

    S.AmU(:,1)    = Au;
    S.PmU(:,:,1)  = Pu;



    for t = 1:nobs
        %       t
        % A = A_t|t-1   & P = P_t|t-1
        
        % Predicted (a priori) state estimate
        A   = T*Au; 
        
        % Predicted (a priori) estimate covariance
        P   = T*Pu*T' + Q; 
        P   =  0.5 * (P+P');
        
        Plarge(:,t) = diag(P);
%         if t == 121
%             fprintf('P');
%             disp(diag(P));
%             fprintf('T');
%             disp(diag(T));
%         end

        % handling the missing data
        [y_t,Z_t,R_t,L_t] = MissData(Y(:,t),Z,R);

        if isempty(y_t)
            Au = A;
            Pu = P;

        else
            PZ  = P*Z_t';
            
            % Innovation (or residual) covariance
            iF  = Z_t*PZ + R_t;
            
            % Optimal Kalman gain
            PZF = PZ/iF;
            
            % Innovation or measurement residual
            V   = y_t - Z_t*A;
            
            % Updated (a posteriori) state estimate
            Au  = A  + PZF * V;
            
            % Updated (a posteriori) estimate covariance
            Pu  = P  - PZF * PZ';
            Pu   =  0.5 * (Pu+Pu');

            if t == T
                disp(t)
                disp(sum(abs(Pu)))
                disp(det(iF))
            end
            
            S.loglik = S.loglik + 0.5*(log(det(inv_FP(iF)))  - V'/iF*V);
        end

        S.Am(:,t)   = A;
        S.Pm(:,:,t) = P;

        % Au = A_t|t   & Pu = P_t|t

        S.AmU(:,t+1)    = Au;
        S.PmU(:,:,t+1)  = Pu;
    end % t

    if isempty(y_t)
        S.KZ = zeros(m,m);
    else
        S.KZ = PZF*Z_t;
    end
end