function [y,Mx,Sx]=center(x)

T=size(x,1);

Mx=nanmean(x);
Sx=nanstd(x);

y=(x-ones(T,1)*Mx)./(ones(T,1)*Sx);