function Res = para_const_FP(X,P, lag)

    % Kalman filter with specified paramaters
    % written for
    % "MAXIMUM LIKELIHOOD ESTIMATION OF FACTOR MODELS ON DATA SETS WITH
    % ARBITRARY PATTERN OF MISSING DATA."
    % by Marta Ba?bura and Michele Modugno

    Z_0 = P.Z_0;
    V_0 = P.V_0;
    A = P.A;
    C = P.C;
    Q = P.Q;
    R = P.R;
    Mx = P.Mx;
    Wx = P.Wx;

    %--------------------------------------------------------------------------
    % Preparation of the data
    %--------------------------------------------------------------------------
    
    [T, N] = size(X);
    
    % Standardise x
    xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);
    y    = xNaN';

    Sf = SKF_FP(y,C,R,A,Q,Z_0,V_0);
    Ss = FIS_FP(y,C,R,A,Q,Sf);

    Pf = Sf.PmU(:,:,2:end);
    Ps = Ss.PmT(:,:,2:end);
    
    Plag{1} = Ps;

    for jk = 1:lag
        for jt = size(Plag{1},3):-1:lag+1;
            As = Pf(:,:,jt-jk)*A'*pinv(A*Pf(:,:,jt-jk)*A'+Q); %Marta corrected from As = Pf(:,:,jt-1)*A'*pinv(A*Pf(:,:,jt-1)*A'+Q);
            Plag{jk+1}(:,:,jt) = As*Plag{jk}(:,:,jt);
        end;
    end;

    Res.Plag = Plag;

    Zsmooth = Ss.AmT;
    Vsmooth = Ss.PmT;

    Zsmooth = Zsmooth';
    x_sm    = Zsmooth(2:end,:)*C';
    X_sm    = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);

    %--------------------------------------------------------------------------
    %   Loading the structure with the results
    %--------------------------------------------------------------------------
    Res.P    = Vsmooth;
    Res.X_sm = X_sm;
    Res.F    = Zsmooth(2:end,:); 