function [z_ecdf] = PITtest(PIT,rvec)
temp = PIT;
z = temp(isnan(temp)==0);
P = length(z);
cumcumz = [];
for r = rvec
    cumcumz = [cumcumz (z<r)-r];
end
v = sum(cumcumz,1)/sqrt(P);

Qv = []; Qvabs = [];
for rs = 1:size(rvec,2)
    Qv = [Qv; v(1,rs)'*v(1,rs)];
    Qvabs = [Qvabs; abs(v(1,rs)')];
end

KvSTnaive = max(Qvabs);
CVMvSTnaive = mean(Qv);


for r = 1:size(rvec,2)
    z_ecdf(:,r) = mean(z < rvec(:,r));
end