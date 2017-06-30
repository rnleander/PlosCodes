% This file test the ability of the code to fit the parameters for 
% the one-stage, one-stage with lag and two-stage models using artificial
% data. The absolute error in the mean and standard deviation of each part 
% of the cell cycle are displayed in plots as the number of data points
% increases.  

clear all

twostagefit=0;

twostagelag=0;

onestagelag=0;

onestagefit=0;

threestagefit=1;

%compute the first three moments of the data for generating an initial
%guess for the paramters



%-------------------------------------------------------------------------

%for fitting parameters using the method of moments
%[x0_d,l0_d,res_d,ERR_d]=moments_method_2stage(dmso);
%[x0_e,l0_e,res_e,ERR_e]=moments_method_2stage(erlot);
%[x0_c,l0_c,res_c,ERR_c]=moments_method_2stage(chx);
%[x0_e]=moments_method_2stage([.0788 .0238 .9003 2.2686],erlot);

%------------------------------------------------------------------------

if onestagefit == 1
    
    data_sizes=[200 1000 10000];
    n=length(data_sizes);
    mu=.05;
    sigma=.08;
    
    mu_error=zeros(n,1);
    sigma_error=zeros(n,1);
    
    p=zeros(n,2);
    ld=zeros(n,1);
    
    %number of data sets to fit for each data size
    k=5;
    
    for i=1:n
        num=data_sizes(i);
        pp=zeros(k,2);
        for j=1:k
    
        data=testdata_IVG(data_sizes(i),mu,sigma);
    

        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        %use true parameters as initial guess
        x0=[mu,sigma];
        [pp(j,:),conf1]=mle(data,'pdf',@onestagepdf2,'start',x0, 'upperbound', [Inf Inf],'lowerbound',[0 0],'options',options)
   
        end
        mu_error(i,1)=sum(abs(mu-pp(:,1)))/k;
        sigma_error(i,1)=sum(abs(sigma-pp(:,2)))/k;
        
    end
    
    plot(log10(data_sizes),mu_error./mu,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('realative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu_error_onestage.fig'));
    
    plot(log10(data_sizes),sigma_error./sigma,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma_error_onestage.fig'));

end
%-------------------------------------------------------------------------

%fit data to onestage model with a lag using mle

if onestagelag==1

    data_sizes=[200 1000 10000];
    n=length(data_sizes);
    mu=.1;
    sigma=.3;
    lag=10;
    
    p=zeros(n,3);
    ld=zeros(n,1);
    
    mu_error=zeros(n,1);
    sigma_error=zeros(n,1);
    lag_error=zeros(n,1);
    
    k=5;
    
    for i=1:n
        num=data_sizes(i);
        pp=zeros(k,3);
        
        for j=1:k
            
            data=testdata_IVG_lag(data_sizes(i),mu,sigma,lag);
        
            options = statset('MaxIter',10000, 'MaxFunEvals',10000);

            %use true parameters for initial guess
            x0=[mu,sigma,lag];
            [pp(j,:),conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options)
         
        
        end
        mu_error(i,1)=sum(abs(mu-pp(:,1)))/k;
        sigma_error(i,1)=sum(abs(sigma-pp(:,2)))/k;
        lag_error(i,1)=sum(abs(lag-pp(:,3)))/k;
    end
    
    
    plot(log10(data_sizes),mu_error./mu,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu_error_onestage_lag.fig'));
    
    plot(log10(data_sizes),sigma_error./sigma,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma_error_onestage_lag.fig'));
    
    plot(log10(data_sizes),lag_error./lag,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in mitotic delay','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('lag_error_onestage_lag.fig'));


end
%-------------------------------------------------------------------------

%fit data to twostage model with lag using mle

if twostagelag==1

    data_sizes=[200 1000 10000];
    n=length(data_sizes);
    
    mu1=.08;
    sigma1=.03;
    mu2=.14;
    sigma2=.6;
    lag=1;
    
    
    p=zeros(n,5);
    ld=zeros(n,1);
    
    mu1_error=zeros(n,1);
    sigma1_error=zeros(n,1);
    mu2_error=zeros(n,1);
    sigma2_error=zeros(n,1);
    lag_error=zeros(n,1);
    
    k=5;
    
    for i=1:n
        num=data_sizes(i);
        pp=zeros(k,5);
        
        for j=1:k
            
            data=testdata_IVG_2stagelag(data_sizes(i),mu1,sigma1,mu2,sigma2,lag);
    
            %use true parameters for initial guess
            x0=[mu1 sigma1 mu2 sigma2 lag];
            
        
            options = statset('MaxIter',10000, 'MaxFunEvals',10000);

           f=@(t,m1,s1,m2,s2,l)twostagepdf_lag(t,m1,s1,m2,s2,l,.01,10^(-8));
            [pp(j,:),conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0],'options',options)
            l=twostagepdf_lag(data,pp(j,1),pp(j,2),pp(j,3),pp(j,4),pp(j,5),.01,10^(-8));
            l=sum(log(l));
            
        
        end
        mu1_error(i,1)=sum(abs(mu1-pp(:,1)))/k;
        sigma1_error(i,1)=sum(abs(sigma1-pp(:,2)))/k;
        mu2_error(i,1)=sum(abs(mu2-pp(:,3)))/k;
        sigma2_error(i,1)=sum(abs(sigma2-pp(:,4)))/k;
        lag_error(i,1)=sum(abs(lag-pp(:,5)))/k;
        
    end
    
    
    plot(log10(data_sizes),mu1_error./mu1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu1_error_twostagelag.fig'));
    
    plot(log10(data_sizes),sigma1_error./sigma1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma1_error_twostagelag.fig'));
    
    plot(log10(data_sizes),mu2_error./mu2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu2_error_twostagelag.fig'));
    
    plot(log10(data_sizes),sigma2_error./sigma2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma2_error_twostagelag.fig'));
    
    plot(log10(data_sizes),lag_error./lag,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in lag','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('lag_error_twostagelag.fig'));

end

%-------------------------------------------------------------------------

%-------------------------------------------------------------------------

%fit data to twostage model using mle

if twostagefit==1

    data_sizes=[200 1000 10000];
    n=length(data_sizes);
    
%       DMSO values
%       E=3.7245e-04
%       h=5.0000e-04
%       l= -665.8195

% 
%     mu1=.0788;
%     sigma1=.0239;
%     mu2=.9087;
%     sigma2=2.3140;

%       Erlotnib values
%       E=1.1591e-04
%       h=5.0000e-04
%       l= -815.4143


    mu1=.0819;
    sigma1=.0152;
    mu2=.1765;
    sigma2=.5875;
    
    p=zeros(n,4);
    ld=zeros(n,1);
    
    mu1_error=zeros(n,1);
    sigma1_error=zeros(n,1);
    mu2_error=zeros(n,1);
    sigma2_error=zeros(n,1);
    
    k=10;
    
    for i=1:n
        num=data_sizes(i);
        pp=zeros(k,4);
        
        for j=1:k
            
            data=testdata_IVG_2stage(data_sizes(i),mu1,sigma1,mu2,sigma2);
    
            %use true parameters for initial guess
            x0=[mu1 sigma1 mu2 sigma2];
            
        
            options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');

           f=@(t,m1,s1,m2,s2)convolv_2invG_adapt_nov(t,m1,s1,m2,s2,.01);
            [pp(j,:),conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)
            l=convolv_2invG_adapt_nov(data,pp(j,1),pp(j,2),pp(j,3),pp(j,4),.01);
            l=sum(log(l));
            
        
        end
        mu1_error(i,1)=sum(abs(mu1-pp(:,1)))/k;
        sigma1_error(i,1)=sum(abs(sigma1-pp(:,2)))/k;
        mu2_error(i,1)=sum(abs(mu2-pp(:,3)))/k;
        sigma2_error(i,1)=sum(abs(sigma2-pp(:,4)))/k;
        
    end
    
    
    plot(log10(data_sizes),mu1_error./mu1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu1_error_twostage.fig'));
    
    plot(log10(data_sizes),sigma1_error./sigma1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma1_error_twostage.fig'));
    
    plot(log10(data_sizes),mu2_error./mu2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu2_error_twostage.fig'));
    
    plot(log10(data_sizes),sigma2_error./sigma2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma2_error_twostage.fig'));

end

%-------------------------------------------------------------------------

%fit data to three-stage model using mle

if threestagefit==1

    data_sizes=[200 1000 10000];
    n=length(data_sizes);
    

% test parameters erlotinib
% 
    mu1=.1879;
    sigma1=.6733;
    mu2=.4182;
    sigma2=.2555;
    mu3=.0982;
    sigma3=.0002;
    
%test parameters MCF
%     mu1=.2532;
%     sigma1=1.0166;
%     mu2=.1139;
%     sigma2=.0548;
%     mu3=.1485;
%     sigma3=.0715;

    
    mu1_error=zeros(n,1);
    sigma1_error=zeros(n,1);
    mu2_error=zeros(n,1);
    sigma2_error=zeros(n,1);
    mu3_error=zeros(n,1);
    sigma3_error=zeros(n,1);
    
    k=10;
    flag=zeros(n,k);
    
    
    ld=zeros(n,k);
    R_error=zeros(n,k);
    ld_check=zeros(n,k);
    
    
    for i=1:n
        num=data_sizes(i);
        pp=zeros(k,6);
        
        for j=1:k
            
            data=testdata_IVG_3stage(data_sizes(i),mu1,sigma1,mu2,sigma2,mu3,sigma3);
            max(data)
            min(data)
    
            %use true parameters for initial guess
            x0=[mu1 sigma1 mu2 sigma2 mu3 sigma3];
            
            options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');

           g=@(t,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(t,m1,s1,m2,s2,m3,s3,.01);
            [pp(j,:),conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options)
            [l,flag(i,j)]=convolv_3invG_nov(data,pp(j,1),pp(j,2),pp(j,3),pp(j,4),pp(j,5),pp(j,6),.01);
            ld(i,j)=sum(log(l));
        
        end
        mu1_error(i,1)=sum(abs(mu1-pp(:,1)))/k;
        sigma1_error(i,1)=sum(abs(sigma1-pp(:,2)))/k;
        mu2_error(i,1)=sum(abs(mu2-pp(:,3)))/k;
        sigma2_error(i,1)=sum(abs(sigma2-pp(:,4)))/k;
        mu3_error(i,1)=sum(abs(mu3-pp(:,5)))/k;
        sigma3_error(i,1)=sum(abs(sigma3-pp(:,6)))/k;
        
        
    
end
        
    
    
    plot(log10(data_sizes),mu1_error./mu1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu1_error_threestage.fig'));
    
    plot(log10(data_sizes),sigma1_error./sigma1,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma1_error_threestage.fig'));
    
    plot(log10(data_sizes),mu2_error./mu2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu2_error_threestage.fig'));
    
    plot(log10(data_sizes),sigma2_error./sigma2,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma2_error_threestage.fig'));
    
    plot(log10(data_sizes),mu3_error./mu3,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \mu','FontSize', 16);
    %axis([2,5, min(mu_error), max(mu_error)]);
    saveas(gcf,sprintf('mu3_error_threestage.fig'));
    
    plot(log10(data_sizes),sigma3_error./sigma3,'ro')
    xlabel('data size','FontSize', 16);
    ylabel('relative error in \sigma','FontSize', 16);
    %axis([2,5, min(sigma_error), max(sigma_error)]);
    saveas(gcf,sprintf('sigma3_error_threestage.fig'));

end





