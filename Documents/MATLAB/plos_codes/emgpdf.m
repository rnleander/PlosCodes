function Y=emgpdf(X,l,m,s)


Y=(l/2)*erfc((-X+m+l*s^2)/(s*2^(1/2))).*exp((l/2)*(-2*X+2*m+l*s^2));
Y=max(realmin, Y);
end


