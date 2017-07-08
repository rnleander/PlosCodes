function Y=onestagepdf2(t,mu,s)

Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));
%The pdf may return values that are zero to witin machine error
%these values are also replaced by realmin
Y=max(Y, realmin);
end
