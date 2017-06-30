%this function evaluates a shifted inverse gaussian distribution at X
function Y=onestagepdf_lag(X,m,s,l)
%m=mu
%s=sigma
%l=lag

Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l)));
%If some cpmponent of X is less than l, the pdf will return imaginary values, 
%in this case we stipilate that the pdf returns realmin
Y=real(Y);
Y=max(Y, realmin);
end