function Y=onestagepdf_prime(t,m,s)
% find the value of the first derivitive of the inverse gaussian
% distribution with parameters m and s at each point in t
% and return a list of these values corresponding to each point t

Y=exp((-(1-m*t).^2)./(2*s^2*t));
Y=(-3/2*t.^(-5/2)+(t.^(-2)/(2*s^2)-m^2/(2*s^2)).*t.^(-3/2))*(1/(s*(2*pi)^.5)).*Y;
end
