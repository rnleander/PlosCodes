function M=gp_max(m,s)
%this function finds the maximum value of the first derivative of an
%inverse Gaussian distribution

a=m^4/(4*s^4);
b=3*m^2/(2*s^2);
c=15/4-m^2/(2*s^4);
d=-5/(2*s^2);
e=1/(4*s^4);

z=roots([a b c d e]);
z=z(imag(z)==0);
z=z(z>=0);
M=onestagepdf_prime(z,m,s);
M=max(abs(M));

end
