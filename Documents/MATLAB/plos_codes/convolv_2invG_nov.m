%This function evaluates the convolution of two inverse Gaussian
%distributions at vector t.

%If the variance in one of the distributions is very small so that the 
%distribution is close the a Dirac delta function, the convolution
%is approximated as a shifted inverse Gaussian, that is 
%one part of the cell cycle is treated as deterministic in length.  
%In this case, the shift or lag is the average time to complete 
%the part of the cycle with the smallest sigma.

%t is a vector of times to divide (or times to complete two parts of the
%cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1,
%s2=sigma2.

function [P,flag]=convolv_2invG_nov(t,m1,s1,m2,s2,h)
flag=0;

E=Inf;
n=length(t);
m=[m1 m2];
s=[s1 s2];
m=max(0,m);
s=max(0,s);
v=(s.^2)./(m.^3);
[v,I]=sort(v);
sd=v.^.5;
m=m(I);
s=s(I);

T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));

%When called from convolv_3invG all of the components of t may be negative.
%In this case, the probability of each is set to realmin

if max(t)<=0
    P=realmin.*ones(size(t));
else
    

    Maxt=max(t);
   
    x=0:h:Maxt;
    
    y=onestagepdf2(x,m(1),s(1));
    z=onestagepdf2(x,m(2),s(2));
    
    if sd(1)<.01
    
        %to estimate the error (in calculating probabilities the of the data), 
        %that results from approximating the first pdf as a point-mass
        %distribtion, find the maximum of the absolute value of the 
        %derivative of the second pdf.
        
    gp=gp_max(m(2),s(2));
    
    %determine the radius, r, of a small interval over which the second pdf, g, is
    %approximately constant and so that g(t) is small for t<r.  
    
    r=min(eps/(3*gp),T2);
    checkval=onestagepdf2(r,m(2),s(2));
    while checkval>=eps/2
        r=r/2;
        checkval=onestagepdf2(r,m(2),s(2));
    end
    
    %get the average value of the first pdf.  This is the point at which
    %its mass is concentrated.
    nu=1/m(1);
    
    %get the maximum value f the second pdf.
    gm=onestagepdf2(T2,m(2),s(2));
    %get upper limit of integral for approximating
    %int_{r+nu}^{infty}f(s)ds.
    Tu=max(100,nu+r+1000*sd(1));
    
        
        checkerror=100;
        hh=.001;
        check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1)));
        while checkerror>10^(-4)
            hh=.5*hh;
            ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1)));
            checkerror=abs(check1-ck1);
            check1=ck1;
        end
   check2=gm*check1;
    if  check2<=eps/3
        
        flag=1;
        sigma=s(2);
        mu=m(2);
        l=1/m(1);
    
    P=onestagepdf_lag(t,mu,sigma,l);
    
    else
    
% find the discrete convolution of the vectors y and z
% the (i-1)th element of v approximates the convolution of the pdfs 
% over [.001, x(i)] as a left-hand Riemann sum.
    C=conv(z,y)*h;
    N=length(y);
% only the first N elements of the convolution are valid
    C=C(1:N);
    I=zeros(1,n);
    P=zeros(1,n);
    for i=1:n
%find element of x that is closest to t(i)
        [~,I(i)]=min((t(i)-x).^2);
        %If t(i)<0 the probability is set to zero, otherwise the
        %probability is approxiated as a value from the vector x.
        if t(i)>0 && I(i)>1
            P(i)=C(I(i)-1);
        else
            P(i)=realmin;
        end
            
    end
    %toc
    P=P';
    
    logP=sum(log(P));
    end
    
    else
        
        % find the discrete convolution of the vectors y and z
% the (i-1)th element of v approximates the convolution of the pdfs 
% over [.001, x(i)] as a left-hand Riemann sum.
    C=conv(z,y)*h;
    N=length(y);
% only the first N elements of the convolution are valid
    C=C(1:N);
    I=zeros(1,n);
    P=zeros(1,n);
    for i=1:n
%find element of x that is closest to t(i)
        [~,I(i)]=min((t(i)-x).^2);
        %If t(i)<0 the probability is set to zero, otherwise the
        %probability is approxiated as a value from the vector x.
        if t(i)>0 && I(i)>1
            P(i)=C(I(i)-1);
        else
            P(i)=realmin;
        end
            
    end
    %toc
    P=P';
    
    logP=sum(log(P));
    end
        
        
end
