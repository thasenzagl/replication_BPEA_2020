function [ ySA, trend, seasonality ] = timeVar_trig_sa( y, S )
%TIMEVAR_TRIG_SA removes the seasonality component from y
%
%   This function models in a Kalman Filter / Kalman Smoother routine
%   the seasonality component and trend of y.
%
%   A. The seasonality is estimated with a time varying trigonometric model 
%   B. The trend with a local level trend model
%
%
%   =======================================================================
%   MODEL
%   =======================================================================
%   
%   State space form:
%   
%   y_{t} = B*F_{t}   + eps_{t} ~ N(0, varIdio)
%   F_{t} = A*F_{t-1} + u_{t}   ~ N(0, [varSeasonVec, varTrend] * eye(S+1))
%   
%   Where:
%
%   F_{t}         = [gamma_{t}; gamma_{t}^{+}; mu_{t}]
%   B             = blkdiag(B_{gamma}, B_{mu})
%   A             = blkdiag(A_{gamma}, A_{mu})
%   varSeasonVec  = [varSeason_{gamma, 1}, varSeason_{gamma^{+}, 1}, ...]
%   varIdio       = 1e-04
%
%
%   =======================================================================
%   A. Time varying trigonometric seasonality
%   =======================================================================
%   
%   We obtain the stochastic trigonometric seasonal component gamma_{t}
%   by having the following of equations:
%
%   gamma_{t} = sum_{j=1}^{S/2} (gamma_{j, t})
%
%   And:
%
%   C_{j} = | gamma_{j, t}     | = | cos(lambda_{j})    sin(lambda_{j}) | * 
%           | gamma_{j, t}^{+} |   | -sin(lambda_{j})   cos(lambda_{j}) |  
%
%           * | gamma_{j, t-1}     | + N(0, varSeason * eye(2))
%             | gamma_{j, t-1}^{+} |
%
%   Where:
%
%   lambda_{j} = 2*pi*j/S;
%   gamma_{j, t} and gamma_{j, t}^{+} are treated as unknown coefficients
%
%   Thus:
%
%   B_{gamma} = (1, 0, ..., 1, 0) in order to impose the initial sum
%   A_{gamma} = blkdiag(A_{gamma}, C_{j}) for each j
%   varSeason = 1
%
%
%   See: Messy Time Series: A Unified Approach, Harvey, Koopman and Penzer
%
%
%   =======================================================================
%   B. Local level trend
%   =======================================================================
%
%   The following equation is used to model a random walk trend for y:
%
%   mu_{t} = mu_{t-1} + N(0, varTrend)   
%
%   Thus:
%
%   B_{mu}   = 1
%   A_{mu}   = 1
%   varTrend = 1
%
%
%   =======================================================================
%   INPUT
%   =======================================================================
%
%   y             NSA variable                                        (1xT)
%   S             No. of seasons (i.e. 12 for monthly variables)      (1x1)
%
%  
%   =======================================================================
%   OUTPUT
%   =======================================================================
%   
%   ySA           Seasonally adjusted y                               (1xT)
%   trend         Trend component                                     (1xT)
%   seasonality   Seasonality component                               (1xT)
%
%
%   =======================================================================
%   AUTHOR
%   =======================================================================
%
%   Filippo Pellegrino, Economist at Now-Casting Economics
%                       Email: filippo.pellegrino@now-casting.com


    %% Initial settings

    C0 = [];
    Q0 = [];
    P0 = [];

    varIdios0  = 1e-04;
    varTrend0  = 1; 
    varSeason0 = 1;


    %% A. Time varying trigonometric seasonality

    for j=1:S/2
        l_j = 2*pi*j/S;
        C_j =  [cos(l_j), sin(l_j); ...
               -sin(l_j), cos(l_j)];

        Q_j = eye(2)*varSeason0;                                                % Eye = Serially and mutually indipendent
        P_j = reshape(pinv(eye(4)-kron(C_j,C_j))*Q_j(:), 2, 2);

        % Store the parameters
        C0 = blkdiag(C0, C_j);
        Q0 = blkdiag(Q0, Q_j);
        P0 = blkdiag(P0, P_j);
    end
    clear l_j C_j Q_j P_j;

    A0  = [1, 0];
    A0  = repmat(A0, 1, S/2);
    F0 = zeros(12, 1);


    %% Local level trend

    C0(end+1, end+1) = 1;
    Q0(end+1, end+1) = varTrend0;
    P0(end+1, end+1) = 0;
    A0 = [A0, 1]; 
    F0 = [F0; y(1)];


    %% Idiosyncratic component
    R0 = varIdios0;


    %% KF / KS
    [FmT, PmT, PmT_1, loglik] = runKF(y, C0, A0, Q0, R0, F0, P0);

    seasonality = A0(1, 1:end-1)*FmT(1:end-1, :);
    seasonality = seasonality(:, 2:end);
    trend       = FmT(end, 2:end);
    ySA         = y-seasonality;

end

