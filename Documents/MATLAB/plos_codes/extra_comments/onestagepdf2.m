function Y=onestagepdf2(t,mu,s)
% find the value of the inverse gaussian with parameters mu and s
% at each point in t and return a list Y with the value cooresponding
% to the points in t

Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));

%The pdf may return values that are zero to witin machine error
%these values are also replaced by realmin
Y=max(Y, realmin);

end
