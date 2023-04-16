function [ A, C, Q, R, initZ, initV] = InitCond_VAR(x,r,p,blocks,optNaN,Rcon,q,NQ,i_idio)


pC = size(Rcon,2);
ppC = max(p,pC);
n_b = size(blocks,2);

OPTS.disp=0;

[xBal,indNaN] = remNaNs_spline(x,optNaN);
[T,N] = size(xBal);
NM = N-NQ;

xNaN = xBal;
xNaN(indNaN) = nan;
C = [];
A = [];
Q = [];
initV = [];

res = xBal;
resNaN = xNaN;
indNaN(1:pC-1,:) = true;
zZz=[];
for i = 1:n_b
    r_i = r(i);
    %--------------------------------------------------------------------------
    % Observation equation
    %--------------------------------------------------------------------------
    C_i = zeros(N,r_i*ppC);
    idx_i = find(blocks(:,i));
    idx_iM = idx_i(idx_i<NM+1);
    idx_iQ = idx_i(idx_i>NM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%MOD MICHELE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if size(idx_iM,1)>1
        [ v, d ] = eigs(cov(res(:,idx_iM)),r_i,'lm',OPTS);
        C_i(idx_iM,1:r_i) = v;
        f = res(:,idx_iM)*v;
    else
        
        C_i(idx_iM,1:r_i) = 1;
        f = res(:,idx_iM);
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = [];
    for kk = 0:max(p+1,pC)-1
        F = [F f(pC-kk:end-kk,:)];
    end
    Rcon_i = kron(Rcon,eye(r_i));
    q_i = kron(q,zeros(r_i,1));
    ff = F(:,1:r_i*pC);
    for j = idx_iQ'
        xx_j = resNaN(pC:end,j);
        if sum(~isnan(xx_j)) < size(ff,2)+2
            xx_j = res(pC:end,j);
        end
        ff_j = ff(~isnan(xx_j),:);
        xx_j = xx_j(~isnan(xx_j));
        iff_j = inv(ff_j'*ff_j);
        Cc = iff_j*ff_j'*xx_j;
        Cc = Cc - iff_j*Rcon_i'*inv(Rcon_i*iff_j*Rcon_i')*(Rcon_i*Cc-q_i);
        C_i(j,1:pC*r_i)=Cc';
    end
    ff = [zeros(pC-1,pC*r_i);ff];
    res = res - ff*C_i';
    resNaN = res;
    resNaN(indNaN) = nan;
    C = [C C_i];

    %--------------------------------------------------------------------------
    % Transition equation
    %--------------------------------------------------------------------------
    z = F(:,1:r_i);
    Z = F(:,r_i+1:r_i*(p+1));
    zZz=[zZz z];
    A_i = zeros(r_i*ppC,r_i*ppC)';
    A_temp = inv(Z'*Z)*Z'*z;
    A_i(1:r_i,1:r_i*p) = A_temp';
    A_i(r_i+1:end,1:r_i*(ppC-1)) = eye(r_i*(ppC-1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_i = zeros(ppC*r_i,ppC*r_i);
    e = z  - Z*A_temp;         % VAR residuals
    Q_i(1:r_i,1:r_i) = cov(e); % VAR covariance matrix

    initV_i = reshape(inv(eye((r_i*ppC)^2)-kron(A_i,A_i))*Q_i(:),r_i*ppC,r_i*ppC);

    A = blkdiag(A,A_i);
    Q = blkdiag(Q,Q_i);
    initV = blkdiag(initV,initV_i);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%MOD MICHELE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = zZz;
r = size(zZz,2);
Z = [];
for kk = 1:p
    Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
end;
z = z(p+1:end,:);

A_temp = (inv(Z'*Z)*Z'*z)';
A(pC+1,1:p)=A_temp(2,[1 3]);
A(1,pC+1:pC+p)=A_temp(1,[2 4 ]);
e = z  - Z*A_temp';         % VAR residuals
Q_temp = cov(e);
Q(pC+1,1)=Q_temp(2,1);% qui era l' errore. c' era: Q(pC+1,1:p)=Q_temp(2,1);
Q(1,pC+1)=Q_temp(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = diag(nanvar(resNaN));

eyeN = eye(N);
eyeN(:,~i_idio) = [];
% Initial conditions
C=[C eyeN];

ii_idio = find(i_idio);
n_idio = length(ii_idio);
B = zeros(n_idio);
S = zeros(n_idio);

for i = 1:n_idio;
    R(ii_idio(i),ii_idio(i)) = 1e-04;

    res_i = resNaN(:,ii_idio(i));
    % number of leading zeros
    leadZero = max( find( (1:T)' == cumsum(isnan(res_i)) ) );
    endZero = max( find( (1:T)' == cumsum(isnan(res_i(end:-1:1))) ) );

    res_i = res(:,ii_idio(i));
    res_i(end-endZero:endZero) = [];
    res_i(1:leadZero) = [];

    BM(i,i) = inv(res_i(1:end-1)'*res_i(1:end-1))*res_i(1:end-1)'*res_i(2:end,:);
    SM(i,i) = cov(res_i(2:end)-res_i(1:end-1)*B(i,i));
end
initViM = diag(1./diag(eye(size(BM,1))-BM.^2)).*SM;


C = [C [zeros(NM,5*NQ);kron(eye(NQ),[1 2 3 2 1])]];
Rdiag = diag(R);
sig_e = Rdiag(NM+1:N)/19;
Rdiag(NM+1:N) = 1e-04;
R = diag(Rdiag);
    %%%%%%%%%%%%%%%%%%%%%%%%%%MOD MICHELE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R(9,9)=1e-100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0 = 0.1;
BQ = kron(eye(NQ),[[rho0 zeros(1,4)];[eye(4),zeros(4,1)]]);
temp = zeros(5);
temp(1,1) = 1;
SQ = kron(diag((1-rho0^2)*sig_e),temp);

initViQ = reshape(inv(eye((5*NQ)^2)-kron(BQ,BQ))*SQ(:),5*NQ,5*NQ);

A = blkdiag(A, BM, BQ);
Q = blkdiag(Q, SM, SQ);

% Initial conditions
initZ = zeros(size(A,1),1); 
initV = blkdiag(initV, initViM, initViQ);