%this function evaluates the convolution of two inverse Gaussian
%distributions with a deterministic lag at vector t
function [P,flag]=twostagepdf_lag(t,m1,s1,m2,s2,l,h,tol)

flag=0;


%organize parameters to check for small values of sigma
%assuming that at most one sigma values is small.

m=[m1 m2];
s=[s1 s2];
v=(s.^2)./(m.^3);
[v,I]=sort(v);
sd=v.^.5;
m=m(I);
s=s(I);

x1=max(t);
x=0:h:500;
%When called from convolv_3invG all of the components of t may be negative.
%In this case, the probability of each is set to realmin

if max(t)<=0
    P=realmin;
else
    
%integrate each pdf
y1=sum(h*onestagepdf2(x,m(1),s(1)));
y2=sum(h*onestagepdf2(x,m(2),s(2)));
%determine the squared errors between the integrals and 1.
err1=(y1-1)^2;
err2=(y2-1)^2;

%if the pdfs do not integrate to 1 within the specified tolerance, use the
%shifted inverse Gaussian instead if the convolution


if err1>tol && err2<=tol
    
        flag=1;
        
        sigma=s(2);
        mu=m(2);
        l2=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
        
        P=onestagepdf_lag(t,mu,sigma,l+l2);

elseif err1<=tol && err2>tol
    
        flag=1;
        
        sigma=s(1);
        mu=m(1);
        l2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
        
        P=onestagepdf_lag(t,mu,sigma,l+l2);
%otherwise, perform the convolution

elseif err1>tol && err2>tol
    
        flag=2;
        
        l2=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)))+(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
        P=zeros(1,length(t));
        P(t==l+l2)=1;
else
    [P,flag]=convolv_2invG_small_sigma_test_march(t'-l,m(1),s(1),m(2),s(2),h,tol);
end
P=max(realmin,P);
end
end