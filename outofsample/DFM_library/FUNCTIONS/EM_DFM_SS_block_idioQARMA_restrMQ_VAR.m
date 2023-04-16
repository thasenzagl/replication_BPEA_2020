function Res = EM_DFM_SS_block_idioQARMA_restrMQ_VAR(X,Par,Res_old)

if nargin<3
    % Standardise x
    Mx = nanmean(X);
    Wx = (nanstd(X));
    
else
    Mx = Res_old.Mx;
    Wx = Res_old.Wx;
end

nQ = Par.nQ;

nM = size(X,2)-nQ;

thresh = 1e-3;
r = Par.r;
p = Par.p;
max_iter = Par.max_iter;

i_idio = logical([ones(nM,1);zeros(nQ,1)]);
i_idio(9) = logical(0);


R_mat = [2 -1 0 0 0;...
    3 0 -1 0 0;...
    2 0 0 -1 0;...
    1 0 0 0 -1]; % matrix of constraints on the loading of quarterly data

q = zeros(4,1);


blocks = Par.blocks;

%--------------------------------------------------------------------------
% Preparation of the data
%--------------------------------------------------------------------------
[T,N] = size(X);


xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);

%--------------------------------------------------------------------------
% Initial Conditions
%--------------------------------------------------------------------------

%Removing missing values (for initial estimators)
optNaN.method = 2; % Remove leading and closing zeros
optNaN.k = 3;

if nargin<3
    [A, C, Q, R, Z_0, V_0] = InitCond_VAR(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
else
    A = Res_old.A;
    C = Res_old.C;
    Q = Res_old.Q;
    R = Res_old.R;
    Z_0 = Res_old.Z_0;
    V_0 = Res_old.V_0;
end;
    
    
    
% some auxiliary variables for the iterations
previous_loglik = -inf;
num_iter = 0;
LL = -inf;
converged = 0;

% y for the estimation is WITH missing data
y = xNaN';


%--------------------------------------------------------------------------
%THE EM LOOP
%--------------------------------------------------------------------------

%The model can be written as
%y = C*Z + e;
%Z = A*Z(-1) + v
%where y is NxT, Z is (pr)xT, etc

%remove the leading and ending nans for the estimation
optNaN.method = 3;
y_est = remNaNs_spline(xNaN,optNaN)';

while (num_iter < max_iter) & ~converged
    [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = ...
        EMstep_VAR(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks);
    
    C = C_new;
    R = R_new;
    A = A_new;
    Q = Q_new;

    % Checking convergence
    if num_iter>2
    [converged,decrease(num_iter+1)] = em_converged(loglik, previous_loglik, thresh,1);
    end
    
    LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
end

%final run of the Kalman filter
%----------------------------------------------
Zsmooth = runKF_FP(y, A, C, Q, R, Z_0, V_0)';
x_sm = Zsmooth(2:end,:)*C';


Res.X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);
Res.F = Zsmooth(2:end,:);

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
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