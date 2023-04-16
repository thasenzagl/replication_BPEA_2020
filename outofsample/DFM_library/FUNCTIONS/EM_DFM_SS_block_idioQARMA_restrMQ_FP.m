function Res = EM_DFM_SS_block_idioQARMA_restrMQ_FP(X,P,Res_old)
%% EM_DFM_SS_block_idioQARMA_restrMQ
%
%   This function computes the maximum likihood estimates using the EM
%   algorithm. The inputs to the EM_DFM function are as follows.
%   X:
%   Res_old: empty array
%   P: structure with the following fields
%       P.r:          number of factors
%       P.p:          number of lags in the VAR process
%       P.max_iter:   maximum number of iterarions
%       P.i_idio:     nx1 vector of logical 1 (true) for monthly variables and
%                     logical 0 for quarterly variables
%       P.Rconstr:    4x5 matrix of restrictions
%       P.q:          4x1 vector of zeros
%       P.nQ:         number of variables with quarterly frequency
%       P.Series:     names of the individual series
%       P.blocks:     nx3 matrix. P.blocks. If an entry is 1 the corresponding
%                     series is part of the block global, hard or soft, respectively
%
%   The EM_DFM function returns the following output:
%   Res:
%
%   The algorithm consists of the following 2 steps.
%       E-step: Calculate expectations of the log-likelihood conditional on
%               the data using estimates from the previous iteration
%       M-step: Maximize the expected log-likelihood to obtain new
%               parameters
%   Those steps are repeated until the difference of the likelihood
%   functions from the current and previous iterations are smaller than
%   the threshold.

    %% Input Variables
    %

    thresh = 1e-3; % threshold can be changed to 1e-4
    r = P.r; % number of factors
    p = P.p; % number of lags in the VAR process
    max_iter = P.max_iter; 
    i_idio = P.i_idio; %1 for monthly, 0 for quarterly variables
    Rmat = P.Rconstr; % Mariano Murasawa restriction
    Rvec = P.Rconstr_vec; % vector of Mariano Murasawa restriction
    q = P.q; 
    nQ = P.nQ; 
    blocks = P.blocks; 
    

    %% Standardise the data
    %
    %   Mx:     1xn vector. Every entry contains mean of every variable (column of X) ignoring NaN values
    %   Wx:     1xn vector. Every entry contains std of every variable (column of X) ignoring NaN values
    %   xNaN:   Txn matrix. Every variable is demeaned and has a std of 1.

    T = size(X,1);

    Mx = nanmean(X);
    Wx = (nanstd(X));
    xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);

    %% Initial Conditions

    % Removing missing values (for initial estimators)
    optNaN.method = 2; % Remove leading and closing zeros
    optNaN.k = 3;

    % y for the estimation is WITH missing data
    y = xNaN';

    if isempty(Res_old) % usually we run function with inputs (X_new(1:iQ+3,:), P, []) so Res.old is [] and thus empty

        [A, C, Q, R, Z_0, V_0] = InitCond_FP(xNaN, r, p, blocks, optNaN, ...
                                          Rmat, Rvec, q, nQ, i_idio);

        % some auxiliary variables for the iterations
        previous_loglik = -inf;
        num_iter = 0;
        LL = -inf;
        converged = 0;


        %% THE EM LOOP

        %The model can be written as
        %y = C*Z + e;
        %Z = A*Z(-1) + v
        %where y is NxT, Z is (pr)xT, etc

        %remove the leading and ending nans for the estimation
        optNaN.method = 3;
        y_est = remNaNs_spline(xNaN,optNaN)';

        while (num_iter < max_iter) && ~converged
            [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = ...
                EMstep_FP(y_est, A, C, Q, R, Z_0, V_0, r, p, Rmat, Rvec, ...
                       q, nQ, i_idio, blocks);

            C = C_new;
            R = R_new;
            A = A_new;
            Q = Q_new;

            % Checking convergence
            if num_iter>2
                [converged, decrease(num_iter+1)] = ...
                    em_converged(loglik, previous_loglik, thresh,1);
            end

            LL = [LL loglik];
            previous_loglik = loglik;
            num_iter =  num_iter + 1;
        end
        
        
        %-----------------------------------------------
    else
        A=Res_old.A; C=Res_old.C; Q=Res_old.Q; R=Res_old.R; Z_0=Res_old.Z_0; V_0=Res_old.V_0;

    end
    %----------------------------------------------
    %final run of the Kalman filter to get the conditional expectation of
    %the factor
    Zsmooth = runKF_FP(y, A, C, Q, R, Z_0, V_0)';
    x_sm = Zsmooth(2:end,:)*C';


    Res.X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);
    Res.F = Zsmooth(2:end,:);


    %%   Loading the structure with the results

    Res.C = C;
    Res.R = R;
    Res.A = A;
    Res.Q = Q;
    Res.Mx = Mx;
    Res.Wx = Wx;
    Res.Z_0 = Z_0;
    Res.V_0 = V_0;
    Res.r = r;
    Res.p = p;
end