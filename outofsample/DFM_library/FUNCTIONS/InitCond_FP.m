function [ A, C, Q, R, initZ, initV] = InitCond_FP(x, r, p, blocks,  ...
    optNaN, Rmat, Rvec, q, NQ, ...
    i_idio)

    % INITCOND determines the initial conditions for the EM algorithm
    %
    %   INPUT:
    %
    %       x           Variables standardised to mean 0 std 1 (T x n)
    %       r           Number of factors per block (1 x no. of blocks)
    %       p           No. of lags (1 x 1)
    %       blocks      Blocks setting matrix (n x no. blocks)
    %       optNaN      remNaNs_spline argument (1 x 1 struct)
    %       Rmat        Mariano Murasawa restriction (no. blocks x 5)
    %       q           zeros(1, no. blocks)
    %       NQ          No. quarterly variables
    %       i_idio      1 if monthly, 0 otherwise (n x 1)
    %
    %
    %   Last edit: Thomas Hasenzagl, 26rd April 2016


    %% Initial settings

    pC = size(Rmat, 2);                                            % Always 5 with Mariano Murasawa restriction
    ppC = max(p,pC);                                               % Max(lags, 5)
    n_b = size(blocks, 2);                                         % no. blocks
    OPTS.disp=0;                                                   % Eig argument
    [xBal,indNaN] = remNaNs_spline(x,optNaN);                      % It fills the NaNs in x using optNaN as argument (with the default optNaN.method replace missing values after removing leading and closing zeros)
    [T, N] = size(xBal);                                           % T is not the same as size(x, 1) because xBal has no leading and closing NaNs
    NM = N - NQ;                                                   % No. monthly variables
    xNaN = xBal;
    xNaN(indNaN) = nan;                                            % Same as xBal with NaNs. Very similar to x, but leading and closing NaNs were removed.
    C = [];
    A = [];
    Q = [];
    F = [];
    initV = [];
    res = xBal;                                                    % Initialise res
    resNaN = xNaN;                                                 % Initialise resNaN
    indNaN(1:pC-1,:) = true;                                       % treat the first pC-1 rows as if they were rows of NaNs


    %% Principal Components of Monthly Variables for each block
    for i = 1:n_b  % loop over each block

        r_i = r(i);     % Number of factors per block
        C_i = zeros(N, r_i*ppC);    % Zero matrix
        idx_i  = find(blocks(:, i));
        idx_iM = idx_i(idx_i < NM + 1); % Index of monthly variables
        idx_iQ = idx_i(idx_i > NM);     %Index of quarterly variables


        if size(idx_iM, 1) > 1
            % Calculate the symmetric positive semidefinite covariance matrix of
            % the data
            Sigma = cov(res(:,idx_iM));

            % Obtain the eigenvectors of the covariance matrix
            [ v, ~ ] = eigs(Sigma,r_i,'lm',OPTS);

            % Loadings are equal to the eigenvectors
            C_i(idx_iM,1:r_i) = v;

            % Obtain the principal components: Use the eigenvectors as new basis
            % vectors and project the data onto the new basis.
            f = res(:,idx_iM)*v;
        else
            % If there is only one monthly series, that series is the factor
            C_i(idx_iM,1:r_i) = 100000/99999;
            f = res(:, idx_iM);
        end

        % Lags of the first principal components. Number of lags is max(p, Pc).
        for kk = 0:ppC-1
            F = [F f(pC-kk:end-kk,:)];
        end

        %% Quarterly Variables
        Rcon_i = kron(Rmat, eye(r_i)); % Mariano Murasawa restriction for block i
        q_i = kron(q, zeros(r_i, 1));   % vector of zeros
        %ff = F(:, 1:r_i*pC); % Resize F, given r (if r=1, ff=F)
        ff = F(:, 1+(i-1)*ppC:ppC*(r_i+i-1));

        for j = idx_iQ' % Iterate for each quarterly variable in block i
            xx_j = resNaN(pC:end, j);
            if sum(~isnan(xx_j)) < size(ff,2)+2
                xx_j = res(pC:end, j); % use the quarterly variable without NaNs
            end

            % remove the months without releases from the quarterly variable
            % and the factor
            ff_j = ff(~isnan(xx_j),:);
            xx_j = xx_j(~isnan(xx_j));

            % regression: xx_j = Cc * iff_j + e
            iff_j = ff_j'*ff_j;
            Cc = iff_j\ff_j'*xx_j; % Regression coefficients
            % Apply the Mariano Murasawa restriction on the parameters
            Cc = Cc - iff_j\Rcon_i'/(Rcon_i/iff_j*Rcon_i')*(Rcon_i*Cc-q_i);
            C_i(j, 1:pC*r_i)=Cc';
        end

        ff  = [zeros(pC-1,pC*r_i); ff];
        res = res - ff*C_i'; % Residuals
        resNaN = res;
        resNaN(indNaN) = nan; % Residuals with NaN on missing observations
        C = [C C_i];  % Store the i-th coefficients


        %% State equation: AR(5) model
        % z_t = (A1 A2 0 0 0) (z_t-1 z_t-2 z_t-3 z_t-4 z_t-5)' + e_t
        %       and A_temp = (A1 A2)

        % State-space notation of the AR(2) (write it as first order VAR)
        %   State Equation:
        %       z_t         A1 A2 0 0 0   z_t-1   e_t
        %       z_t-1       1  0  0 0 0   z_t-2    0
        %       z_t-2   =   0  1  0 0 0 * z_t-3 +  0
        %       z_t-3       0  0  1 0 0   z_t-4    0
        %       Z_t-4       0  0  0 1 0   z_t-5    0
        %   Observation equation: z_t = (1 0 0 0 0)*alpha_t
        %       where alpha_t = (z_t z_t-1 z_t-2 z_t-3 z_t-4)'

        z = F(:, 1+(i-1)*ppC:r_i+(i-1)*ppC);
        Z = F(:, r_i+1+(i-1)*ppC:r_i*(p+1)+(i-1)*ppC);

        A_i = zeros(r_i*ppC,r_i*ppC)';
        A_temp = inv_FP(Z'*Z)*Z'*z; % estimate A_temp = (A1 A2) with OLS
        A_i(1:r_i,1:r_i*p)= A_temp'; % Initialize A_i
        A_i(r_i+1:end,1:r_i*(ppC-1)) = eye(r_i*(ppC-1));

        Q_i = zeros(ppC*r_i,ppC*r_i); % Initialise Q_i
        e = z  - Z*A_temp; % VAR residuals
        Q_i(1:r_i,1:r_i) = cov(e); % VAR covariance matrix

        % covariance matrix of state vector z
        % vec[cov(z)] = (I_ppc^2 - kron(A,A))^(-1) vec[cov(e)]
        % reshape vec[cov(z)] to get covariance matrix
        initV_i = reshape(inv_FP(eye((r_i*ppC)^2)-kron(A_i,A_i))*Q_i(:),r_i*ppC,r_i*ppC);

        % Store the ith parameters of the transition equation
        A = blkdiag(A,A_i);
        Q = blkdiag(Q,Q_i);
        initV = blkdiag(initV,initV_i); % blockdiagonal

    end

    %% Idiosyncratic Component

    R = diag(nanvar(resNaN));
    eyeN = eye(N);
    eyeN(:,~i_idio) = [];

    C=[C eyeN];

    ii_idio = find(i_idio);
    n_idio = length(ii_idio);
    BM = zeros(n_idio);
    SM = zeros(n_idio);

    % - Monthly variables (idiosyncratic component)
    for i = 1:n_idio;
        R(ii_idio(i), ii_idio(i)) = 1e-04;

        res_i = resNaN(:,ii_idio(i));

        % Number of leading zeros
        leadZero = max( find( (1:T)' == cumsum(isnan(res_i)) ) );
        endZero  = max( find( (1:T)' == cumsum(isnan(res_i(end:-1:1))) ) );

        res_i                  = res(:, ii_idio(i));
        %res_i(end-endZero:endZero) = [];
        res_i(end-endZero+1:end) = [];
        res_i(1:leadZero)      = [];

        BM(i,i) = (res_i(1:end-1)' * res_i(1:end-1))\res_i(1:end-1)'*res_i(2:end,:); % AR(1) coefficient for the idiosyncratic component (monthly variables)
        
        if BM(i,i) >= 1
            BM(i,i) = eps();
        end
        
        SM(i,i) = cov(res_i(2:end) - res_i(1:end-1)*BM(i,i)); % Cov matrix for the observation equation, taking in account the AR(1) idiosyncratic component (monthly variables)
    end

    initViM = diag(1./diag(eye(size(BM,1))-BM.^2)).*SM;

    C = [C [zeros(NM,size(Rmat, 2)*NQ); kron(eye(NQ), Rvec)]];
    Rdiag = diag(R);

    % - Quarterly variables (idiosyncratic component)
    sig_e = Rdiag(NM+1:N)/19; %/19;
    Rdiag(NM+1:N) = 1e-04;
    R = diag(Rdiag);


    rho0 = 0.1;
    BQ = kron(eye(NQ),[[rho0 zeros(1, size(Rmat, 1))]; [eye(size(Rmat, 1)), zeros(size(Rmat, 1),1)]]);             % AR(1) coefficient for the idiosyncratic component (quarterly variables)
    temp = zeros(size(Rmat, 2));
    temp(1,1) = 1;
    SQ = kron(diag((1-rho0^2)*sig_e), temp); % Cov matrix for the observation equation, taking in account the AR(1) idiosyncratic component (quarterly variables)

    initViQ = reshape(inv_FP(eye((size(Rmat, 2)*NQ)^2)-kron(BQ,BQ))*SQ(:),size(Rmat, 2)*NQ,size(Rmat, 2)*NQ);

    A = blkdiag(A, BM, BQ);
    Q = blkdiag(Q, SM, SQ);

    % Initial conditions
    initZ = zeros(size(A,1),1);
    initV = blkdiag(initV, initViM, initViQ);
end