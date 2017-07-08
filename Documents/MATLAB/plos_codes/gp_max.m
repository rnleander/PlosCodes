function M=gp_max(m,s)
%this function finds the maximum value of the first derivative of an
%inverse Gaussian distribution with parameters m and s

% the second derivitive of inverse gaussian can be represented as a
% polynomial with the following coefficents, a through e

a=m^4/(4*s^4);
b=3*m^2/(2*s^2);
c=15/4-m^2/(2*s^4);
d=-5/(2*s^2);
e=1/(4*s^4);

% find the roots of the second derivitive
z=roots([a b c d e]);

% we are only interested in the positive real roots
%z=z(imag(z)==0);
z=z(abs(imag(z))<=0.000000000000001);
z=z(z>=0);

% find the first derivitive at each of these roots
M=onestagepdf_prime(z,m,s);

% return the maximu value of the first derivative
M=max(abs(M));

end
