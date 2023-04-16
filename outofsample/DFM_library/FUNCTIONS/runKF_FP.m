function [AmT, PmT, PmT_1, loglik] = runKF_FP(Y, T, Z, Q, R, A_0, P_0)
    
    % RUNKF Run Kalman Filter and Smoother
    %
    %   MODEL:
    %
    %       Y_t   = Z * A_t + epsilon_t
    %       A_t+1 = T * A_t + u_t
    %
    %
    %   INPUT:
    %
    %       Y          Data                                      (nobs x n)
    %       T          Estimated coefficients                    (nobs x n)
    %       Z          Factor loadings                           (nobs x n)
    %       Q          u_t ~ N(0, Q)                             (nobs x n)
    %       R          epsilon_t ~ N(0, R)                       (nobs x n)
    %       A_0        Initial value of the state vector A       (nobs x n)
    %       P_0        Initial covariance of the state vector A  (nobs x n)
    %
    %
    %   OUTPUT:
    %
    %       AmT        Estimated A_t|T                           (nobs x n)
    %       PmT        P_t|T = Cov(A_t|T)                        (nobs x n)
    %       PmT_1      P_{t-1}|T = Cov(A_{t-1}|T)                (nobs x n)
    %       loglik     u_t ~ N(0, Q)                             (nobs x n)
    %
    %       Note: m is the dim of state vector and t = 1 ...T is time
    %
    %
    %   Last edit: Filippo Pellegrino, Now-casting Economics Ltd, 6/11/15
    
    S = SKF_FP(Y, Z, R, T, Q, A_0, P_0);
    S = FIS_FP(Y, Z, R, T, Q, S);

    AmT = S.AmT;
    PmT = S.PmT;
    PmT_1 = S.PmT_1;
    loglik = S.loglik;