function [L,varargout] = kalcvf(data, lead, a, F, b, H, var, varargin)
%KALCVF The Kalman filter
%
%   State space model is defined as follows:
%     z(t+1) = a+F*z(t)+eta(t)     (state or transition equation)
%       y(t) = b+H*z(t)+eps(t)     (observation or measurement equation)
%
%   [logl, <pred, vpred, <filt, vfilt>>]=kalcvf(data, lead, a, F, b, H, var, <z0, vz0>)
%   computes the one-step prediction and the filtered estimate, as well as their covariance matrices.
%   The function uses forward recursions, and you can also use it to obtain k-step estimates.
%
%   The inputs to the KALCVF function are as follows:
%     data is a Ny×T matrix containing data (y(1), ... , y(T)).
%     lead is the number of steps to forecast after the end of the data.
%        a is an Nz×1 vector for a time-invariant input vector in the transition equation.
%        F is an Nz×Nz matrix for a time-invariant transition matrix in the transition equation.
%        b is an Ny×1 vector for a time-invariant input vector in the measurement equation.
%        H is an Ny×Nz matrix for a time-invariant measurement matrix in the measurement equation.
%      var is an (Ny+Nz)×(Ny+Nz) matrix for a time-invariant variance matrix for
%             the error in the transition equation and the error in the measurement equation,
%             that is, [eta(t)', eps(t)']'.
%       z0 is an optional Nz×1 initial state vector.
%      vz0 is an optional Nz×Nz covariance matrix of an initial state vector.
%
%   The KALCVF function returns the following output:
%     logl is a value of the average log likelihood function of the SSM
%             under assumption that observation noise eps(t) is normally distributed
%     pred is an optional Nz×(T+lead) matrix containing one-step predicted state vectors.
%    vpred is an optional Nz×Nz×(T+lead) matrix containing mean square errors of predicted state vectors.
%     filt is an optional Nz×T matrix containing filtered state vectors.
%    vfilt is an optional Nz×Nz×T matrix containing mean square errors of filtered state vectors.
%
%   The initial state vector and its covariance matrix of the time invariant Kalman filters
%   are computed under the stationarity condition:
%          z0 = (I-F)\a
%         vz0 = (I-kron(F,F))\V(:)
%   where F and V are the time invariant transition matrix and the covariance matrix of transition equation noise,
%   and vec(V) is an Nz^2×1 column vector that is constructed by the stacking Nz columns of matrix V.
%   Note that all eigenvalues of the matrix F are inside the unit circle when the SSM is stationary.
%   When the preceding formula cannot be applied, the initial state vector estimate is set to a
%   and its covariance matrix is given by 1E6I. Optionally, you can specify initial values.
%
%   This is a M-file for MATLAB.
%   Copyright 2002-2003 Federal Reserve Bank of Atlanta
%   $Revision: 1.2 $  $Date: 2003/03/19 19:16:17 $
%   Iskander Karibzhanov 5-28-02.
%   Master of Science in Computational Finance
%   Georgia Institute of Technology
%==========================================================================
% Revision history:
%
%  03/19/2003  -  algorithm and interface were adapted from SAS/IML KALCVF subroutine for use in MATLAB M file
%
%==========================================================================

T = size(data,2);
Nz = size(a,1);
Ny = size(b,1);

nin = nargin;
if nin~=7 && nin~=9
   error('Seven or nine input arguments required.')
end
if nin==9
   z = varargin{1};
   P = varargin{2};
end

nout = nargout;
if nout~=1 && nout ~=3 && nout ~=5
   error('One, three, or five output arguments required.')
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
if nin==9 && any(size(z)~=[Nz 1])
   error('z0 must be column vector of length Nz')
end
if nin==9 && any(size(P)~=[Nz Nz])
   error('vz0 must be Nz by Nz matrix')
end

% V(t) and R(t) are variances of eta(t) and eps(t), respectively,
% and G(t) is a covariance of eta(t) and eps(t)
V = var(1:Nz,1:Nz);
R = var(Nz+1:end,Nz+1:end);
G = var(1:Nz,Nz+1:end);

if nin==7
   e = eig(F);
   if all(all(e*e'-eye(Nz)))
      z = (eye(Nz)-F)\a;
      P = reshape((eye(Nz^2)-kron(F,F))\V(:),Nz,Nz);
   else
      z = a;
      P = eye(Nz)*1e6;
   end
end
if nout>1
   pred = zeros(Nz,T+lead);
   pred(:,1) = z;
   vpred = zeros(Nz,Nz,T+lead);
   vpred(:,:,1) = P;
end
if nout==5
   filt = zeros(Nz,T);
   vfilt = zeros(Nz,Nz,T);
end

L = 0;

for t=1:T
   D = H*P*H'+R;
   dy = data(:,t)-H*z-b;
   ddy = D\dy;
   L = L+log(det(D))+dy'*ddy;
   if nout==5
      PH = P*H';
      filt(:,t) = z+PH*ddy;
      vfilt(:,:,t) = P-PH/D*PH';
   end
   if t<T || lead>0
      FP = F*P;
      FPHG = FP*H'+G;
      z = F*z+FPHG*ddy+a;
      P = FP*F'-FPHG/D*FPHG'+V;
      if nout>1
         pred(:,t+1) = z;
         vpred(:,:,t+1) = P;
      end
   end
end
if lead>1 && nout>1
   for t=T+2:T+lead
      z = F*z+a;
      P = F*P*F'+V;
      pred(:,t) = z;
      vpred(:,:,t) = P;
   end
end

L = -(Ny*log(2*pi)+L/T)/2;

if nout>1
   varargout(1) = {pred};
   varargout(2) = {vpred};
end
if nout==5
   varargout(3) = {filt};
   varargout(4) = {vfilt};
end