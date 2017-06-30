% This function evaluates the convolution of two inverse Gaussian
% distributions at vector t.

% If the variance in one of the distributions is very small so that the 
% distribution is close the a Dirac delta function, the convolution
% is approximated as a shifted inverse Gaussian, that is, 
% one part of the cell cycle is treated as deterministic in length.  
% In this case, the shift or lag is the average time to complete 
% the part of the cycle that has the smallest standard deviation.

% t is a vector of times to divide (or times to complete two parts of the
% cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1,
% s2=sigma2.

function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
% Input parameters:
% t = the list of points at which to evaluate the distribution
% m1 = mu for the first distribution
% s1 = sigma for the first distribution
% m2 = mu for the second distribution
% s2 = sigma for the second distribution
% h = step size

% Outputs:
% P = probability at each point corresponding to points in t
% h = final step size
% flag = indicates if we used the Dirac delta
% E = is the relative error in the likelihood of the data due to the numerical integration

% log the parameters we were called with
%fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ... 
%    m1, s1, m2, s2);

% ????
flag=0;

% a constant that is used in determining if the Dirac approximation should be applied.
eps=.01;

% ????
E=Inf;

% number of points to be evaluated
n=length(t);

% store m and s for each sub-distribution in a list
m=[m1 m2];
s=[s1 s2];

% make sure we dont have negative values for m or s by replacement with
% zero if we do
m=max(m,0);
s=max(s,0);

% find the variance for both sub-distributions
v=(s.^2)./(m.^3);

% find the standard deviation for each sub-distribution
sd=v.^.5;

% reorder m and s so that the sub-distribution with the smallest
% variance comes first.  So the fist part might be approximated as a Dirac delta.
[v,I]=sort(v);
m=m(I);
s=s(I);


% T1 appears to be unused, should probably remove
% T2 is only used inside FUNCTION CHECK_APPROXIMATABLE
% and should probably be moved inside that scope
% T2 is the mode of the second distribution, determined by a precomputed analytic expression.
T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));

%When called from convolv_3invG all of the components of t may be negative.
%In this case, the probability of each is set to realmin
if max(t)<=0
    P=realmin.*ones(size(t));
    
% otherwise we need to calculate P for each point in t
else
    % find the largest point in t
    Maxt=max(t);
   
    % produce a range of evenly spaced points to evaluate at between
    % 0 and Maxt with step size h. the even spacing is important when
    % calculating the convolution later
    x=0:h:Maxt;
    
    % evaluate the sub-distributions at each point in x
    y=onestagepdf2(x,m(1),s(1));
    z=onestagepdf2(x,m(2),s(2));
    
    % if the first pdf is very concentrated, check to see if it can be
    % approximated as a point-mass distribution
    if sd(1)<.01
    
        % BEGIN FUNCTION TailMass
        % Input parameters: m, s, eps, T2, nu, sd,
        % Outputs: check2
            % to estimate the error (in calculating probabilities the of the
            % data), that results from approximating the first pdf as a
            % point-mass distribtion, find the maximum of the absolute value
            % of the derivative of the second pdf.
            gp=gp_max(m(2),s(2));

            % determine the radius, r, of a small interval over which the
            % 1. second pdf, g, is approximately constant, i.e. changes by less than eps/3
            % over any interval with that radius
            % and 2. g(t) is small for t<r.  ?????
            r=min(eps/(3*gp),T2);
            checkval=onestagepdf2(r,m(2),s(2));
            while checkval>=eps/2
                r=r/2;
                checkval=onestagepdf2(r,m(2),s(2));
            end

            % get the average value of the first pdf.  This is the point at
            % which its mass is concentrated.
            nu=1/m(1);

            % get the maximum value of the second pdf,f.
            gm=onestagepdf2(T2,m(2),s(2));

            % ????
            % Tu is the upper limit of integral for approximating
            % int_{r+nu}^{infty}g(s)ds. This is in latex.
            Tu=max(100,nu+r+1000*sd(1));

            % ????
            checkerror=100;
            hh=.001;
            numinusr=nu-r;
            nuplusr=nu+r;
            teeyou=Tu;
            %we are integrating the first pdf from zero to nu-r, and from nu+r to infinity to see if the tails
            %are small in probability.
            LeftTail=.001*sum(onestagepdf2((0:.001:numinusr),m(1),s(1)));
            RightTail=.001*sum(onestagepdf2((nuplusr:.001:teeyou),m(1),s(1)));
            check1=LeftTail+RightTail;
            %Reduce step size in above Riemann sum until the error is small,
            %meaning that the Riemann sum is converging.
            while checkerror>10^(-4)
                hh=.5*hh;
                ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1)));
                checkerror=abs(check1-ck1);
                check1=ck1;
            end
            check2=gm*check1;
            %     end
        % END FUNCTION CHECK_APPROXIMATABLE
        
        if  check2<=eps/3
            %If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1.
            %Flag is set to 1 to indicate this.
            flag=1;
            sigma=s(2);
            mu=m(2);
            %l for lag.
            l=1/m(1);
            P=onestagepdf_lag(t,mu,sigma,l);
        else
            % BEGIN FUNCTION DOTHECONVOLUTION_OUTER
            % Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x
            % Outputs: P
            
                % BEGIN FUNCTION ApproxConvolv
                % Input parameters: z, y, h, n, t, i, I, x
                % Outputs: logP0
                    % find the discrete convolution of the vectors y and z
                    % the (i-1)th element of v approximates the convolution of the pdfs 
                    % over [0, x(i)] as a left-hand Riemann sum.
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
                    P0=max(realmin,P);
                    logP0=sum(log(P0));
                % END FUNCTION DOTHECONVOLUTION_INNER
            %Keep reducing the step size in the numerical integration until we are happy with the error.
                while E>=.001*abs(logP0)
                    h1=.5*h;
                    x=0:h1:Maxt;
                    y=onestagepdf2(x,m(1),s(1));
                    z=onestagepdf2(x,m(2),s(2));

                    % BEGIN FUNCTION DOTHECONVOLUTION_INNER
                    % Input parameters: z, y, h1, n, t, i, I, x
                    % Outputs: logP1
                        % find the discrete convolution of the vectors y and z
                        % the (i-1)th element of v approximates the convolution of the pdfs 
                        % over [.001, x(i)] as a left-hand Riemann sum.
                        C=conv(z,y)*h1;
                        N=length(y);
                        % only the first N elements of the convolution are valid
                        C=C(1:N);
                        I=zeros(1,n);
                        P=zeros(1,n);
                        for i=1:n
                            %find element of x that is closest to t(i)
                            [~,I(i)]=min((t(i)-x).^2);
                            %If t(i)<0 the probability is set to zero, otherwise the
                            %probability is approximated as a value from the vector x.
                            if t(i)>0 && I(i)>1
                                P(i)=C(I(i)-1);
                            else
                                P(i)=realmin;
                            end
                        end
                        %toc
                        P=P';
                        P1=max(realmin,P);
                        logP1=sum(log(P1));
                    % END FUNCTION DOTHECONVOLUTION_INNER

                    E=abs(logP1-logP0);
                    P0=P1;
                    logP0=logP1;
                    h=h1;
                end
                P=P0;
            % END FUNCTION DOTHECONVOLUTION_OUTER
        end   
    % pdf is not very concentrated so compute the convolution directly
    else
        % BEGIN FUNCTION DOTHECONVOLUTION_OUTER
        % Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x
        % Outputs: P
            % BEGIN FUNCTION DOTHECONVOLUTION_INNER
            % Input parameters: z, y, h, n, t, i, I, x
            % Outputs: logP0
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
                P0=max(realmin,P);
                logP0=sum(log(P0));
            % END FUNCTION DOTHECONVOLUTION_INNER

            while E>=.001*abs(logP0)
                h1=.5*h;
                x=0:h1:Maxt;
                y=onestagepdf2(x,m(1),s(1));
                z=onestagepdf2(x,m(2),s(2));

                % BEGIN FUNCTION DOTHECONVOLUTION_INNER
                % Input parameters: z, y, h1, n, t, i, I, x
                % Outputs: logP1
                    % find the discrete convolution of the vectors y and z
                    % the (i-1)th element of v approximates the convolution of the pdfs 
                    % over [.001, x(i)] as a left-hand Riemann sum.
                    C=conv(z,y)*h1;
                    N=length(y);
                    % only the first N elements of the convolution are valid
                    C=C(1:N);
                    I=zeros(1,n);
                    P=zeros(1,n);
                    for i=1:n
                    %find element of x that is closest to t(i)
                        [~,I(i)]=min((t(i)-x).^2);
                        %If t(i)<0 the probability is set to zero, otherwise the
                        %probability is approximated as a value from the vector x.
                        if t(i)>0 && I(i)>1
                            P(i)=C(I(i)-1);
                        else
                            P(i)=realmin;
                        end
                    end
                    %toc
                    P=P';
                    P1=max(realmin,P);
                    logP1=sum(log(P1));
                % END FUNCTION DOTHECONVOLUTION_INNER

                E=abs(logP1-logP0);
                P0=P1;
                logP0=logP1;
                h=h1;
            end
            P=P0;
        % BEGIN FUNCTION DOTHECONVOLUTION_OUTER
    end
end
end
