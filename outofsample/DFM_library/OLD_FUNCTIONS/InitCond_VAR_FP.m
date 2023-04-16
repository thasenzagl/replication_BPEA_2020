function [ A, C, Q, R, initZ, initV] = InitCond_FPv2(x, r, p, blocks, ...
                                       optNaN, Rmat, Rvec, q, NQ, i_idio)
   
    % INITCOND determines the initial conditions for the EM algorithm
    %
    %   INPUT:
    %       
    %       x           Variables standardised to mean 0 std 1 (T x n)
    %       r           Number of factors per block (1 x no. of blocks)
    %       p           No. of lags (1 x 1)
    %       blocks      Blocks setting matrix (n x no. blocks)
    %       optNaN      remNaNs_spline argument (1 x 1 struct)
    %       Rcon        Mariano Murasawa restriction (no. blocks x 5)
    %       q           zeros(1, no. blocks)
    %       NQ          No. quarterly variables
    %       i_idio      1 if monthly, 0 otherwise (n x 1)
    %
    %
    %   Last edit: Filippo Pellegrino, 23rd December 2015
    
    
    %% Initial settings
    
    pC = size(Rmat, 2);                                                     % Always 5 with Mariano Murasawa restriction
    ppC = max(p,pC);                                                        % Max(lags, 5)
    n_b = size(blocks, 2);                                                  % no. blocks   
    
    OPTS.disp=0;                                                            % Eig argument

    [xBal,indNaN] = remNaNs_spline(x,optNaN);                               % It fills the NaNs in x using using optNaN as argument (with the default optNaN removes the leading and the closing NaN rows)
    [T, N] = size(xBal);                                                    % T is not the same as size(x, 1) because xBal has not leading and closing NaNs             
    NM = N - NQ;                                                            % No. monthly variables

    xNaN = xBal;                                                            
    xNaN(indNaN) = nan;                                                     % Same as xBal with NaNs
    
    C = [];
    A = [];
    Q = [];
    initV = [];

    res = xBal;                                                             % Initialise res
    resNaN = xNaN;                                                          % Initialise resNaN
    indNaN(1:pC-1,:) = true;                                                
    
    
    %% Iterate for each block
    
    z_mat = [];
    Z_mat = [];
    
    for i = 1:n_b
        
        r_i = r(i);                                                         % No. of factors for block i
        
        
        % =================================================================
        % Observation equation
        % =================================================================
        
        C_i = zeros(N, r_i*ppC);                                            % zeros(n, no. factors for block i times max(lags, 5))
        idx_i  = find(blocks(:, i));                                        % general index: 1 if a variable belong to block i, 0 otherwise
        idx_iM = idx_i(idx_i < NM + 1);                                     % monthly index: 1 if a monthly variable belong to block i, 0 otherwise
        idx_iQ = idx_i(idx_i > NM);                                         % quarterly index: 1 if a quarterly variable belong to block i, 0 otherwise
        
        
        % -> Monthly variables
        [ v, ~ ] = eigs(cov(res(:,idx_iM)),r_i,'lm',OPTS);                  % v is the loading matrix - aka: a matrix whose columns are the eigenvectors corresponding to the largest r eigenvalues (for the monthly variables of block i) 
        C_i(idx_iM,1:r_i) = v;                                              
        f = res(:,idx_iM)*v;                                                % f are the factors calculated using PCA (see also Maximum Likelihood Estimation of Factor Models on Data Sets with Arbitrary Pattern of Missing Data, page 9, note 8)
        F = [];
        for kk = 0:max(p+1,pC)-1                                            % F is f lagged with max(p, Pc) lags
            F = [F f(pC-kk:end-kk,:)];
        end
        
        % -> Quarterly variables
        Rcon_i = kron(Rmat, eye(r_i));                                      % Mariano Murasawa restriction for block i (given the associated no. of factors) 
        q_i = kron(q, zeros(r_i, 1));                                       % By default zeros(no. blocks, 1)
        ff = F(:, 1:r_i*pC);                                                % Resize F, given r (note: if r = 1, ff = F)
        
        for j = idx_iQ'                                                     % Iterate for each quarterly variable which belong to the i block
            xx_j = resNaN(pC:end, j);                                       % It starts from the first non-NaN quarterly value
            if sum(~isnan(xx_j)) < size(ff,2)+2 
                xx_j = res(pC:end, j);                                      % Adjust xx_j if sum(~isnan(xx_j)) < size(ff,2)+2 
            end
            ff_j = ff(~isnan(xx_j),:);                                      % It skips the NaNs, in order to take the quarterly value of the factors accordingly
            xx_j = xx_j(~isnan(xx_j));                                      % It does the same with variable j
            iff_j = inv(ff_j'*ff_j);
            Cc = iff_j*ff_j'*xx_j;                                          % Regression coefficients for the j-th variable, in the observational equation
            Cc = Cc - iff_j*Rcon_i'*inv(Rcon_i*iff_j*Rcon_i')*(Rcon_i*Cc-q_i); % It applies the Mariano Murasawa restriction on the parameters
            
            % =============================================================
            % NOTE - BY FP
            % -> In the "Nowcasting" appendix the formula for Cc is
            %    different, but the above is correct
            % =============================================================
            
            C_i(j,1:pC*r_i)=Cc';                                            
        end
        
        ff  = [zeros(pC-1,pC*r_i); ff];
        res = res - ff*C_i';                                                % Residuals
        resNaN = res;
        resNaN(indNaN) = nan;                                               % Residuals with NaN on missing observations
        C = [C C_i];                                                        % Store the i-th coefficients

        % =================================================================
        % Transition equation (FP MOD)
        % =================================================================
        
        % Store z and Z
        
        z = F(:,1:r_i);                                                     % It takes the factors at time t
        Z = F(:,r_i+1:r_i*(p+1));                                           % It takes the p lags for the factors
        
        z_mat = [z_mat, z];
        Z_mat = [Z_mat, Z];
    end

    A_tmp = inv(Z_mat'*Z_mat)*Z_mat'*z_mat;
    A_tmp = A_tmp';                                                         % => size = size(z, 1) * size(Z, 2)                                                       
    z_hat = Z_mat*A_tmp';
    e     = cov(z_mat - z_hat);
    
    A     = zeros(sum(r, 2)*ppC, sum(r, 2)*ppC);
    Q     = zeros(sum(r, 2)*ppC);
    
    % -> indexes
    endInd_A       = cumsum(r*ppC, 2);
    startInd_A     = endInd_A - r*ppC + 1;
    adjEndInd_A    = endInd_A - r*p;
    rp_per_block   = endInd_A - adjEndInd_A;     
    endInd_A       = startInd_A + rp_per_block - 1;
    range_A        = sort([startInd_A(:); endInd_A(:)]);
    
    % -> move from A_tmp to A and from Q_tmp to Q
    A(startInd_A(:), range_A)       = A_tmp;
    Q(startInd_A(:), startInd_A(:)) = e;
    
    for i=1:n_b
        r_i  = r(i);                                                         % No. of factors for block i        
        sizeEyeM = r_i*(ppC-1);
        eyeM     = eye(sizeEyeM);
        A(endInd_A(i):endInd_A(i)+sizeEyeM-1, ...
          startInd_A(i):startInd_A(i)+sizeEyeM-1) = eyeM;                
    end
    
    initV = reshape( inv(eye((size(A, 1))^2)-kron(A, A))*Q(:), ...
                     size(A, 1), size(A, 1));
    

    %% Last estimation
         
    % =====================================================================
    % Observation equation
    % =====================================================================
    
    R = diag(nanvar(resNaN));                                               
    eyeN = eye(N);                                          
    eyeN(:,~i_idio) = [];
    
    C=[C eyeN];                                                                                                                          
    
    ii_idio = find(i_idio);
    n_idio = length(ii_idio);
    B = zeros(n_idio);
    S = zeros(n_idio);

    % - Monthly variables (idiosyncratic component)
    for i = 1:n_idio;                                       
        R(ii_idio(i), ii_idio(i)) = 1e-04;

        res_i = resNaN(:,ii_idio(i));
        
        % Number of leading zeros
        leadZero = max( find( (1:T)' == cumsum(isnan(res_i)) ) );
        endZero  = max( find( (1:T)' == cumsum(isnan(res_i(end:-1:1))) ) );

        res_i                      = res(:, ii_idio(i));
        res_i(end-endZero:endZero) = [];
        res_i(1:leadZero)          = [];

        BM(i,i) = inv(res_i(1:end-1)'*res_i(1:end-1))*res_i(1:end-1)'*res_i(2:end,:); % AR(1) coefficient for the idiosyncratic component (monthly variables)
        SM(i,i) = cov(res_i(2:end)-res_i(1:end-1)*B(i,i));                            % Cov matrix for the observation equation, taking in account the AR(1) idiosyncratic component (monthly variables)
    end
    initViM = diag(1./diag(eye(size(BM,1))-BM.^2)).*SM;

    C = [C [zeros(NM,size(Rmat, 2)*NQ); kron(eye(NQ), Rvec)]];
    Rdiag = diag(R);
    
    
    % =====================================================================
    % NOTE - BY FP
    % - > Why 19?
    % =====================================================================
    
    % - Quarterly variables (idiosyncratic component)
    sig_e = Rdiag(NM+1:N)/19;                                               
    Rdiag(NM+1:N) = 1e-04;
    R = diag(Rdiag);
    rho0 = 0.1;
    BQ = kron(eye(NQ),[[rho0 zeros(1, size(Rmat, 1))]; [eye(size(Rmat, 1)), zeros(size(Rmat, 1),1)]]);             % AR(1) coefficient for the idiosyncratic component (quarterly variables)
    temp = zeros(size(Rmat, 2));
    temp(1,1) = 1;
    SQ = kron(diag((1-rho0^2)*sig_e), temp);                                 % Cov matrix for the observation equation, taking in account the AR(1) idiosyncratic component (quarterly variables)

    initViQ = reshape(inv(eye((size(Rmat, 2)*NQ)^2)-kron(BQ,BQ))*SQ(:),size(Rmat, 2)*NQ,size(Rmat, 2)*NQ);

    A = blkdiag(A, BM, BQ);
    Q = blkdiag(Q, SM, SQ);

    % Initial conditions
    initZ = zeros(size(A,1),1);
    initV = blkdiag(initV, initViM, initViQ);
end
    