%this function evaluates the convolution of three inverse gaussian
%distributions at vector t
function [P,h,flag,E]=convolv_3invG_nov(t,m1,s1,m2,s2,m3,s3,h)
flag=0;
n=length(t);
Maxt=max(t);
x=0:h:Maxt;
E=Inf;
eps=.01;

m=[m1 m2 m3];
s=[s1 s2 s3];
m=max(0,m);
s=max(0,s);
v=(s.^2)./(m.^3);
[v,I]=sort(v);
sd=v.^.5;
m=m(I);
s=s(I);

T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
T3=(1/m(3))*((1+(9/4)*(s(3)^4/m(3)^2))^.5-(3/2)*(s(3)^2/m(3)));

        

if sd(1)<.01
    
        %to estimate the error (in calculating probabilities the of the data), 
        %that results from approximating the first pdf as a point-mass
        %distribtion, find the maximum of the absolute value of the 
        %derivative of the second pdf.
        
    gp=min(gp_max(m(2),s(2)),gp_max(m(3),s(3)));
   
    %define the radius, r, of a small interval over which the second pdf is
    %nearly constant and close to zero for t<r. 
    T=max(T2,T3);
    r=min(eps/(3*gp),T);
    
    checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3)));
    while checkval>=eps/2
        r=r/2;
        checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3))) ;
    end
    
    %estimate the maximum value of the convolution of the second two pdfs
    gm=min(onestagepdf2(T2,m(2),s(2)), onestagepdf2(T3,m(3),s(3)));
    
    
    %get the average value of the first pdf.  This is the point at which
    %its mass is concentrated.
    nu=1/m(1);
    
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
    if check2<=eps/3
        
        flag=1;
        l=1/m(1);
    
    P=convolv_2invG_adapt_nov(t-l,m(2),s(2),m(3),s(3),h);
    
    else
    
    [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h);
    z=onestagepdf2(x',m(1),s(1));
    
% approximate the convolution of the two pdfs
% the ith element of v gives approximates convolution over [0, x(i+1)]
% as a left-hand Riemann sum

C=conv(z,y)*h;
N=length(y);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(1,n);
P=zeros(1,n);
for i=1:n
    %find element of x that is closest to t(i)
    [M,I(i)]=min((t(i)-x).^2);
    if I(i)>1
        P(i)=C(I(i)-1);
    end
    if I(i)==1
        P(i)=0;
    end
end
P=P';
P0=max(realmin,P);
logP0=sum(log(P0));
    
    while E>=.001*abs(logP0)

        h1=.5*h;
    
        x=0:h1:Maxt;
        
        [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h1);
        z=onestagepdf2(x',m(1),s(1));
    
        % approximate the convolution of the two pdfs
        % the ith element of v gives approximates convolution over [.001, x(i+1)]
        % as a left-hand Riemann sum
        C=conv(z,y)*h1;
        N=length(y);
        % only the first N elements of the convolution are valid
        C=C(1:N);
        I=zeros(1,n);
        P=zeros(1,n);
        for i=1:n
        %find element of x that is closest to t(i)
        [~,I(i)]=min((t(i)-x).^2);
        if I(i)>1
            P(i)=C(I(i)-1);
        end
        if I(i)==1
            P(i)=0;
        end
        end
        P=P';
        P1=max(realmin,P);
        logP1=sum(log(P1));
        
        E=abs(logP1-logP0);
    
        P0=P1;
        logP0=logP1;
        h=h1;

    end
    P=P0;
    logP=sum(log(P));
    end
%otherwise, perform the convolution

else
    flag=0;
    [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h);
    z=onestagepdf2(x',m(1),s(1));
    
% approximate the convolution of the two pdfs
% the ith element of v gives approximates convolution over [.001, x(i+1)]
% as a left-hand Riemann sum
C=conv(z,y)*h;
N=length(y);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(1,n);
P=zeros(1,n);
for i=1:n
    %find element of x that is closest to t(i)
    [M,I(i)]=min((t(i)-x).^2);
    if I(i)>1
        P(i)=C(I(i)-1);
    end
    if I(i)==1
        P(i)=0;
    end
end
P=P';
P0=max(realmin,P);
logP0=sum(log(P0));
    
    while E>=.001*abs(logP0)

        h1=.5*h;
    
        x=0:h1:Maxt;
        
        [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h1);
        z=onestagepdf2(x',m(1),s(1));
    
        % approximate the convolution of the two pdfs
        % the ith element of v gives approximates convolution over [.001, x(i+1)]
        % as a left-hand Riemann sum
        C=conv(z,y)*h1;
        N=length(y);
        % only the first N elements of the convolution are valid
        C=C(1:N);
        I=zeros(1,n);
        P=zeros(1,n);
        for i=1:n
        %find element of x that is closest to t(i)
        [~,I(i)]=min((t(i)-x).^2);
        if I(i)>1
            P(i)=C(I(i)-1);
        end
        if I(i)==1
            P(i)=0;
        end
        end
        P=P';
        P1=max(realmin,P);
        logP1=sum(log(P1));
        
        E=abs(logP1-logP0);
    
        P0=P1;
        logP0=logP1;
        h=h1;

    end
    P=P0;
end
end