function Y=onestagepdf_prime(t,m,s)
Y=exp((-(1-m*t).^2)./(2*s^2*t));
Y=(-3/2*t.^(-5/2)+(t.^(-2)/(2*s^2)-m^2/(2*s^2)).*t.^(-3/2))*(1/(s*(2*pi)^.5)).*Y;
end
