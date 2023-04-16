function [sm, vsm] = kalcvs(data, a, F, b, H, var, pred, vpred, varargin)
%KALCVS The Kalman smoothing
%
%   State space model is defined as follows:
%     z(t+1) = a+F*z(t)+eta(t)     (state or transition equation)
%       y(t) = b+H*z(t)+eps(t)     (observation or measurement equation)
%
%   [sm, vsm] = kalcvs(data, a, F, b, H, var, pred, vpred, <un, vun>)
%   uses backward recursions to compute the smoothed estimate z(t|T) and its covariance matrix, P(t|T),
%   where T is the number of observations in the complete data set.
%
%   The inputs to the KALCVS function are as follows:
%     data is a Ny×T matrix containing data (y(1), ... , y(T))'.
%        a is an Nz×1 vector for a time-invariant input vector in the transition equation.
%        F is an Nz×Nz matrix for a time-invariant transition matrix in the transition equation.
%        b is an Ny×1 vector for a time-invariant input vector in the measurement equation.
%        H is an Ny×Nz matrix for a time-invariant measurement matrix in the measurement equation.
%      var is an (Ny+Nz)×(Ny+Nz) covariance matrix for the errors in the transition equation and
%                measurement equation noises, that is, [eta(t)', eps(t)']'.
%     pred is an Nz×T matrix containing one-step forecasts (z(1|0), ... , z(T|T-1))'.
%    vpred is an Nz×Nz×T matrix containing mean square error matrices of predicted state vectors (P(1|0), ... , P(T|T-1))'.
%       un is an optional Nz×1 vector containing u(T).
%      vun is an optional Nz×Nz covariance matrix containing U(T).
%
%   The KALCVS function returns the following output:
%       sm is an Nz×T matrix containing smoothed state vectors (z(1|T), ... , z(T|T))'.
%      vsm is an Nz×Nz×T matrix containing covariance matrices of smoothed state vectors (P(1|T), ... , P(T|T))'.
%
%   This is a M-file for MATLAB.
%   Copyright 2002-2003 Federal Reserve Bank of Atlanta
%   $Revision: 1.1 $  $Date: 2003/03/31 21:38:05 $
%   Iskander Karibzhanov 3-31-03.
%   Master of Science in Computational Finance
%   Georgia Institute of Technology
%===============================================================================================================
% Revision history:
%
%  03/27/2003  -  algorithm and interface were adapted from SAS/IML KALCVS subroutine for use in MATLAB M file
%
%===============================================================================================================*/

T = size(data,2);
Nz = size(a,1);
Ny = size(b,1);

nin = nargin;
if nin~=8 && nin~=10
   error('Eight or ten input arguments required.')
end
if nin==10
   u = varargin{1};
   vu = varargin{2};
else
   u = zeros(Nz,1);
   vu = zeros(Nz);
end

nout = nargout;
if nout~=2
   error('Two output arguments required.')
end

% check input matrix dimensions
if size(data,1)~=Ny
   error('data and b must have the same number of rows')
end
if size(a,2)~=1
   error('a must be column vector')
end
if any(size(F)~=[Nz Nz])
   error('F must be square')
end
if size(b,2)~=1
   error('b must be column vector')
end
if any(size(H)~=[Ny Nz])
   error('H must be Ny by Nz matrix')
end
if any(size(var)~=[(Ny+Nz) (Ny+Nz)])
   error('var must be (Ny+Nz) by (Ny+Nz) matrix')
end
if any(size(pred)~=[Nz T])
   error('pred must be Nz by T matrix')
end
if any(size(vpred)~=[Nz Nz T])
   error('vpred must be Nz by Nz by T matrix')
end
if nin==10 && any(size(un)~=[Nz 1])
   error('un must be column vector of length Nz')
end
if nin==10 && any(size(vun)~=[Nz Nz])
   error('vun must be Nz by Nz matrix')
end

% V(t) and R(t) are variances of eta(t) and eps(t), respectively,
% and G(t) is a covariance of eta(t) and eps(t)
V = var(1:Nz,1:Nz);
R = var(Nz+1:end,Nz+1:end);
G = var(1:Nz,Nz+1:end);

sm = zeros(Nz,T);
vsm = zeros(Nz,Nz,T);

for t=T:-1:1
    z = pred(:,t);
    P = vpred(:,:,t);
    D = H*P*H'+R;
    dy = data(:,t)-H*z-b;
    ddy = D\dy;
    K = (F*P*H'+G)/D;
    L = F-K*H;
    u = H'*ddy+L'*u;
    vu = H'*(D\H)+L'*vu*L;
    sm(:,t) = z+P*u;
    vsm(:,:,t) = P-P*vu*P;
end