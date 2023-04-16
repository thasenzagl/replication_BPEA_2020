function S = FIS_FP(Y,Z,R,T,Q,S)
%______________________________________________________________________
% Fixed intervall smoother (see Harvey, 1989, p. 154)
% FIS returns the smoothed state vector AmT and its covar matrix PmT
% Use this in conjnuction with function SKF
%______________________________________________________________________
% INPUT
%        Y         Data                                 (nobs x n)
%        S Estimates from Kalman filter SKF
%          S.Am   : Estimates     a_t|t-1                  (nobs x m)
%          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
%          S.AmU  : Estimates     a_t|t                    (nobs x m)
%          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)
% OUTPUT
%        S Smoothed estimates added to above
%          S.AmT  : Estimates     a_t|T                    (nobs x m)
%          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
%          S.PmT_1 : Cov(a_ta_t-1|T)
%        where m is the dim of state vector and t = 1 ...T is time

    [m nobs] = size(S.Am);


    %% FP, 6th of November: All the squeeze should be removed - They are not doing anything

    % Initialize S.AmT and S.PmT
    S.AmT             = zeros(m,nobs+1);
    S.PmT             = zeros(m,m,nobs+1);
    S.AmT(:,nobs+1)   = squeeze(S.AmU(:,nobs+1));
    S.PmT(:,:,nobs+1) = squeeze(S.PmU(:,:,nobs+1));
    S.PmT_1(:,:,nobs) = (eye(m)-S.KZ) *T*squeeze(S.PmU(:,:,nobs));

    J_2 = squeeze(S.PmU(:,:,nobs)) * T' * pinv(squeeze(S.Pm(:,:,nobs)));

    for t = nobs:-1:1
        PmU = squeeze(S.PmU(:,:,t));
        Pm1 = squeeze(S.Pm(:,:,t));
        P_T = squeeze(S.PmT(:,:,t+1));
        P_T1 = squeeze(S.PmT_1(:,:,t));

        J_1 = J_2;

        S.AmT(:,t)   = S.AmU(:,t) + J_1 * (S.AmT(:,t+1) - T * S.AmU(:,t));
        S.PmT(:,:,t) = PmU        + J_1 * (P_T - Pm1) * J_1';

        if t>1
            J_2 = squeeze(S.PmU(:,:,t-1)) * T' * pinv(squeeze(S.Pm(:,:,t-1)));
            S.PmT_1(:,:,t-1) = PmU*J_2'+J_1*(P_T1-T*PmU)*J_2';
        end
    end
end