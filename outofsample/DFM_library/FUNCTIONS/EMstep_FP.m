function  [C_new, R_new, A_new, Q_new, Z_0_new, V_0_new, loglik] = EMstep_FP(y, A, ...
           C, Q, R, Z_0, V_0, r, p, Rmat, Rvec, q, nQ, i_idio, blocks)


    % EMSTEP executes the EM algorithm
    %
    %   INPUT:
    %       
    %       y           Variables standardised to mean 0 std 1 (T x n)
    %       A           Transition equation coefficients
    %       C           Observation equation coefficients
    %       Q           Transition equation var-cov matrix
    %       R           Observation equation var-cov matrix
    %       Z_0         
    %       V_0         
    %       r           No. factors
    %       p           No. lags
    %       R_mat       Mariano Murasawa restriction
    %       q           
    %       nQ          No. quarterly variables
    %       i_idio      
    %       blocks      
    %
    %
    %   Last edit: Filippo Pellegrino, 28th December 2015

    
    %% Initial settings
    
    [n,T] = size(y);
    nM = n-nQ;
    pC = size(Rmat,2);
    ppC = max(p,pC);
    n_b = size(blocks,2);
    
    
    %% E-Step
    % Compute the (expected) sufficient statistics for a single Kalman filter sequence.
    % Running the Kalman filter with the current estimates of the parameters
    % See BGR "Nowcasting" for more details (Appendix C, pp. 37, penultimate paragraph)
    
    [Zsmooth, Vsmooth, VVsmooth, loglik] = runKF_FP(y, A, C, Q, R, Z_0, V_0);

    A_new   = A;
    Q_new   = Q;
    V_0_new = V_0;

    
    %% M-step
    %  The new parameters theta (j+1) are estimated through the max. of the
    %  expected log-likelihood (from the previous iteration) with respect
    %  to theta (j)
    
    
    %% -> transition equation (Iterate for each block)
   
    for i = 1:n_b
        r_i = r(i); % No. factors ith block                                                        
        rp  = r_i*p; % // times lags                                                        
        rp1 = sum(r(1:i-1))*ppC; % Index: starting point for the factors (taking in account the dimension of the restriction)                                           
        
        vecInd1 = rp1+1:rp1+r_i*ppC;
        vecInd2 = rp1+1:rp1+rp;
        
        A_i = A(vecInd1,vecInd1);% Transition equation coefficients, ith block                                           
        Q_i = Q(vecInd1,vecInd1);% Transition equation var-cov matrix, ith block                                           
        
        Z_t   = Zsmooth(vecInd2, 2:end);
        Z_t_1 = Zsmooth(vecInd2, 1:end-1);
        
        V_t   = sum(Vsmooth(vecInd2,vecInd2, 2:end), 3);
        V_t_1 = sum(Vsmooth(vecInd2,vecInd2, 1:end-1), 3);
        VV_t  = sum(VVsmooth(vecInd2,vecInd2, :), 3);
        
        EZZ    = Z_t*Z_t'     + V_t; % E(Z'Z)                                        
        EZZ_BB = Z_t_1*Z_t_1' + V_t_1; % E(Z(-1)'Z(-1))                                     
        EZZ_FB = Z_t*Z_t_1'   + VV_t; % E(Z'Z(-1))                                      

        A_i(1:r_i,1:rp)  = EZZ_FB(1:r_i,1:rp) / (EZZ_BB(1:rp,1:rp));
        Q_i(1:r_i,1:r_i) = (EZZ(1:r_i,1:r_i) - ...
                            A_i(1:r_i,1:rp)*EZZ_FB(1:r_i,1:rp)') / T;
        
        % Store the expected value of the parameters: 
        A_new(vecInd1,vecInd1) = A_i;
        Q_new(vecInd1,vecInd1) = Q_i;
        V_0_new(vecInd1,vecInd1) = Vsmooth(vecInd1,vecInd1,1);
    end

    
    %% -> idiosyncratic
    
    rp1 = sum(r)*ppC;
    niM = sum(i_idio(1:nM));
    
    Z_eps_t   = Zsmooth(rp1+1:end,2:end);
    Z_eps_t_1 = Zsmooth(rp1+1:end,1:end-1);
    
    V_eps_t   = sum(Vsmooth(rp1+1:end,rp1+1:end, 2:end),3);
    V_eps_t_1 = sum(Vsmooth(rp1+1:end,rp1+1:end, 1:end-1), 3);
    VV_eps_t  = sum(VVsmooth(rp1+1:end,rp1+1:end,:),3);
    
    EZZ       = diag(diag(Z_eps_t * Z_eps_t'))     + diag(diag(V_eps_t));   % E(epsilon'epsilon)
    EZZ_BB    = diag(diag(Z_eps_t_1 * Z_eps_t_1')) + diag(diag(V_eps_t_1)); % E(epsilon(-1)'epsilon(-1))
    EZZ_FB    = diag(diag(Z_eps_t * Z_eps_t_1'))   + diag(diag(VV_eps_t));  % E(epsilon'epsilon(-1))

    A_i = EZZ_FB * diag(1./diag((EZZ_BB)));
    Q_i = (EZZ - A_i*EZZ_FB') / T;
    
    % Store the parameters (idiosyncratic component):
    A_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = A_i(1:niM,1:niM);                  % Transition eq. coefficients             
    Q_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = Q_i(1:niM,1:niM);                  % Transition eq. var-cov matrix
    V_0_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = diag(diag(Vsmooth(rp1+1:rp1+niM,rp1+1:rp1+niM,1)));
    
    
    %% -> observation equation

    Z_0_new = Zsmooth(:,1); 
    nanY = isnan(y);
    y(nanY) = 0;

    C_new = C;                                                              % Observation eq. coefficients

    % =====================================================================
    % Set indexes 
    % =====================================================================
 
    bl = unique(blocks,'rows');                                             % Unique block combinations
    n_bl = size(bl,1);      
    bl_idxM = [];
    bl_idxQ = [];
    R_con = [];
    q_con = [];
    
    % Iterate for each block:
    for i = 1:n_b
        bl_idxQ = [bl_idxQ repmat(bl(:,i),1,r(i)*ppC)];
        bl_idxM = [bl_idxM repmat(bl(:,i),1,r(i)) zeros(n_bl,r(i)*(ppC-1))];
        R_con = blkdiag(R_con, kron(Rmat,eye(r(i))));
        q_con = [q_con;zeros(r(i)*size(Rmat,1),1)];
    end

    bl_idxM = logical(bl_idxM);
    bl_idxQ = logical(bl_idxQ);

    % Idio:
    i_idio_M = i_idio(1:nM);
    n_idio_M = length(find(i_idio_M));
    c_i_idio = cumsum(i_idio);
    
    
    % Iterate for each block:
    for i = 1:n_bl
        bl_i = bl(i,:);
        rs = sum(r(logical(bl_i)));
        idx_i = find(ismember(blocks,bl_i,'rows'));

        % =================================================================
        % Monthly variables
        % =================================================================
        
        idx_iM = idx_i(idx_i<nM+1);
        n_i = length(idx_iM);

        denom = zeros(n_i*rs,n_i*rs);
        nom = zeros(n_i,rs);

        i_idio_i = i_idio_M(idx_iM);
        i_idio_ii = c_i_idio(idx_iM);
        i_idio_ii = i_idio_ii(i_idio_i);
        for t=1:T
            nanYt = diag(~nanY(idx_iM,t));
            denom = denom + kron(Zsmooth(bl_idxM(i,:),t+1)*Zsmooth(bl_idxM(i,:),t+1)'...
                    + Vsmooth(bl_idxM(i,:),bl_idxM(i,:),t+1),nanYt);
            nom = nom + y(idx_iM,t)*Zsmooth(bl_idxM(i,:),t+1)'...%here's the modification
                  - nanYt(:,i_idio_i)*(Zsmooth(rp1+i_idio_ii,t+1)*Zsmooth(bl_idxM(i,:),t+1)'...
                  + Vsmooth(rp1+i_idio_ii,bl_idxM(i,:),t+1));
        end
        vec_C = (denom)\nom(:);
        C_new(idx_iM,bl_idxM(i,:)) = reshape(vec_C,n_i,rs);

        
        % =================================================================
        % Quarterly variables
        % =================================================================

        idx_iQ = idx_i(idx_i>nM);
        rps = rs*ppC;

        R_con_i = R_con(:,bl_idxQ(i,:));
        q_con_i = q_con;
        no_c = ~(any(R_con_i,2));
        R_con_i(no_c,:) = [];
        q_con_i(no_c,:) = [];

        for j = idx_iQ'
            denom = zeros(rps,rps);
            nom = zeros(1,rps);
            idx_jQ = j-nM;
            i_idio_jQ = (rp1+n_idio_M+size(Rmat, 2)*(idx_jQ-1)+1:rp1+n_idio_M+size(Rmat, 2)*idx_jQ);
            V_0_new(i_idio_jQ,i_idio_jQ) = Vsmooth(i_idio_jQ,i_idio_jQ,1);
            A_new(i_idio_jQ(1),i_idio_jQ(1)) = A_i(i_idio_jQ(1)-rp1,i_idio_jQ(1)-rp1);
            Q_new(i_idio_jQ(1),i_idio_jQ(1)) = Q_i(i_idio_jQ(1)-rp1,i_idio_jQ(1)-rp1);

            for t=1:T
                nanYt = diag(~nanY(j,t));
                denom = denom + kron(Zsmooth(bl_idxQ(i,:),t+1)*Zsmooth(bl_idxQ(i,:),t+1)'...
                    +Vsmooth(bl_idxQ(i,:),bl_idxQ(i,:),t+1),nanYt);
                nom = nom + y(j,t)*Zsmooth(bl_idxQ(i,:),t+1)';
                nom = nom -...
                    nanYt*(Rvec*Zsmooth(i_idio_jQ,t+1)*Zsmooth(bl_idxQ(i,:),t+1)'+...
                    Rvec*Vsmooth(i_idio_jQ,bl_idxQ(i,:),t+1));
            end

            C_i = denom\nom';
            C_i_constr = C_i - denom\R_con_i'/(R_con_i/denom*R_con_i')*(R_con_i*C_i-q_con_i);
            
            % =============================================================
            % NOTE - BY FP
            % -> In the "Nowcasting" appendix the formula for C_i_constr is
            %    different, but the above is correct
            % =============================================================

            C_new(j,bl_idxQ(i,:)) = C_i_constr;
        end
    end

    % =====================================================================
    % Observation equation var-cov matrix
    % =====================================================================

    R_new = zeros(n,n);
    
    for t=1:T
        nanYt = diag(~nanY(:,t));
        R_new = R_new + (y(:,t)-nanYt*C_new*Zsmooth(:,t+1))*(y(:,t)-nanYt*C_new*Zsmooth(:,t+1))'...
            +nanYt*C_new*Vsmooth(:,:,t+1)*C_new'*nanYt...
            +(eye(n)-nanYt)*R*(eye(n)-nanYt);
    end

    R_new = R_new/T;
    RR = diag(R_new);
    RR(i_idio_M) = 1e-04;
    RR(nM+1:end) = 1e-04;
    R_new = diag(RR);

end