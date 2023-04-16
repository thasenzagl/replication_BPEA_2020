%--------------------------------------------------------------------------
% KALMAN FILTER
%--------------------------------------------------------------------------
function [xsmooth, Vsmooth] = runKF_lag(y, A, C, Q, R, x_0, Sig_0, k)

    [n,r] = size(C);
    % if k>4
    %     k=3
    % end

    if k>0
        C = single([C zeros(n,k*r)]);
        A = single(blkdiag(A,zeros(k*r)));
        A(r+1:end,1:k*r) = single(eye(k*r));
        Q = single(blkdiag(Q,zeros(k*r)));
        x_0 = single([x_0;zeros(k*r,1)]);
        Sig_0 = single(blkdiag(Sig_0,zeros(k*r,k*r)));
    end


    S = SKF(y,C,R,A,Q,x_0, Sig_0);

    S = FIS(y,C,R,A,Q,S);


    xsmooth = S.AmT(1:r,:);
    Vsmooth= S.PmT;

end