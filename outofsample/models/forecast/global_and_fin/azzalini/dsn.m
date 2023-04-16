function p=dsn(X,mu,sl,sr);
%density of a split normal distribution
A = (sqrt(2*pi)*(sl+sr)/2)^(-1);

left  = exp(-1/2*(X-mu).^2*1/(sl^2));
right = exp(-1/2*(X-mu).^2*1/(sr^2));

p = A*(left.*(X<mu)+right.*(X>=mu));