function imt()
% This file fits EMG, one-, two-, and three-stage stochasttic models to IMT data.

clear all

%This section of code gets the IMT data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for erlotinib
%load('erlot_imts_April2017.mat')
%data=imt_b;

%for AT1 
%load('AT1_imts_April2017.mat')
%data=imt_b;

%for MCF
%load('MCF_imts_April2017.mat')
%data=imt_b;

%for DMSO 
%load('DMSO_imts_April2017.mat')
%data=imt_b;

%for CHX
%load('CHX_imts_April2017.mat')
%data=imt_b;

%for FUCCI data
%GetProcessedDataParts
load('FUCCI_April2017.mat')
%data=imt_b;
data=G2Time_b;
%data=G1Time_b;

%for PC9 cells
%load('PC9_April2017.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose model to fit
twostagefit=1;
onestagelag=0;
onestagefit=2;
threestagefit=0;
emgfit=0;
twostagelag=0;

%get sample statistics for fitting initializing the model parameters
num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));

%Fit the EMG model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if emgfit == 1
    % BEGIN FUNCTION FIT_EMG
    % Input parameters: C1, C2, data
    % Outputs: 
    
        % prepare statistical variables
        vry = [.25 .5 .75]';  
        c1=C1*vry;
        c2=C2*vry;
        %we vary the parameters so the the Gaussian and exponential parts of
        %the cell cycle are responsible for a fraction of the total mean and
        %variance in the IMT.
        lam_v=1./c1;
        mu_v=c1;
        sig_v=c2.^.5;
        N = length(vry);
        
        % prepare parameter seeds
        pp = cell(N^3);
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    pp{N^2*(i-1)+N*(j-1)+k} = [lam_v(i), sig_v(j), mu_v(k)];
                end
            end
        end
        P = zeros(N^3,3);
        for ii = 1:N^3
            P(ii,:) = pp{ii};
        end
        
        % optimize parameters
        ep = zeros(N^3,3);
        le = -realmax*ones(N^3,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        for i=1:length(P)
            x0=P(i,:);
            options = statset('MaxIter',10000, 'MaxFunEvals',10000);
            [ep(i,:),econf]=mle(data,'pdf',@emgpdf,'start',x0, 'lowerbound',[0 0 0],'upperbound',[100 100 100],'options',options);
            econfint(:,:,i)=econf(:);
            l=emgpdf(data,ep(i,1),ep(i,2),ep(i,3));
            le(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_le,ind_le]=max(le);
        ep_max=ep(ind_le,:);
        confint_max=econfint(:,:,ind_ld);
    % END FUNCTION FIT_EMG
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit the one-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onestagefit == 1
    % BEGIN FUNCTION FIT_ONESTAGE
    
        % prepare statistical variables
        mu=1/C1;
        sigma=(C2/C1^3)^.5;
        vry=.5:.5:2;
        m=mu*vry;
        s=sigma*vry;
        N=length(vry);
        
        % maybe these should be moved down into optimize parameters
        % section?
        pd=zeros(N^2,2);
        ld=-realmax*ones(N^2,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        
        % prepare parameter seeds
        pp = cell(N^2);
        for i = 1:N
            for j = 1:N
                pp{(i-1)*N+j} = [m(i), s(j)];
            end
        end
        P = zeros(N^2,2);
        for ii = 1:N^2
            P(ii,:) = pp{ii};
        end
        
        % optimize parameters
        for i=1:N^2
            x0 = P(i,:);
            [p,conf1]=mle(data,'pdf',@onestagepdf2,'start',x0, 'upperbound', [Inf Inf],'lowerbound',[0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            l=onestagepdf2(data,p(1),p(2));
            ld(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        confint_max=confint(:,:,ind_ld);
    % END FUNCTION FIT_ONESTAGE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit one-stage model with lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onestagelag == 1
    % BEGIN FUNCTION FIT_ONESTAGELAG
    
        % prepare statistical parameters
        mu = C3/(3*C2^2);
        sigma = (C3^3/(27*C2^5))^.5;
        lag = C1-3*C2^2/C3;
        vryv=[0.5 1 2];
        vrym=[.25 .5 .75];
        m = mu*vrym;
        s = sigma*vryv;
        lag = lag*vrym;
        N = length(vrym);
        
        % prepare parameter seeds
        pp = cell(N^3);
        for i = 1:N
            for j = 1:N
                for k = 1:N
                pp{N^2*(i-1)+N*(j-1)+k} = [m(i), s(j), lag(k)];
                end
            end
        end
        P = zeros(N^3,3);
        for ii = 1:N^3
            P(ii,:) = pp{ii};
        end
        
        % optimize parameters
        pd = zeros(N^3,3);
        ld = -realmax*ones(N^3,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        for i=1:length(P)
            x0=P(i,:);
            [p,conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            l=onestagepdf_lag(data,p(1),p(2),p(3));
            ld(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        confint_max=confint(:,:,ind_ld);
    % END FUNCTION FIT_ONESTAGELAG
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit two-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if twostagefit == 1
    % BEGIN FUNCTION FIT_TWOSTAGE
    
        % prepare statistical parameters
        vry = [.25 .5 .75]';
        % vrys=[.01 1 10]';
        % [x0, l_moments, min_res] = moments_method_2stage(data);
        % m1 = x0(1)*vry; 
        % s1 = x0(2)*vrys;
        % m2 = x0(3)*vry;
        % s2 = x0(4)*vrys;
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        N = length(vry);

        % prepare parameter seeds
        
        %get all pairs of the form [m(i),s(j)]
        %these pairs represent all possible unique 
        %parameter choices for each part of the cell
        %cycle.  
        pcomb = allcomb(m,s);
        %place paramter pairs into a cell.  The parameters choices for each part
        %are now indexed
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %get all pairs of indices for the parameter 
        %choices for each part of the cycle to get all 
        %parameter choices for the entire cycle
        id = allcomb(1:length(pcomb),1:length(pcomb));
        %sort the pairs in ascending order.  
        %This equates choices of the form [i,j] and [j,i].
        id = sort(id,2);
        %remove repeats
        id = unique(id,'rows');
        %create a matrix of unique parameter choices for the cell cycle
        P = zeros(length(id),4);
        for ii = 1:length(id)
            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}];
        end

        % optimize parameters
        pd=zeros(length(P),4);
        ld = NaN*ones(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        flag=zeros(length(id),1);
        for i=1:length(id)  
            x0 = P(i,:);
            f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
            %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt2(x,m1,s1,m2,s2,.01,4);
            [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            %[l,flag(i)]=convolv_2invG_small_sigma_test_var(data,p(1),p(2),p(3),p(4),.01,4);
            [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(data,p(1),p(2),p(3),p(4),.01);
            l=sum(log(l));
            ld(i)=l    

        end
        
        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        ld_true=zeros(length(ld),1);
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_adapt_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_TWOSTAGE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit two-stage model with lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if twostagelag == 1
    % BEGIN FUNCTION FIT_TWOSTAGELAG
    
        % prepare statistical parameters
        vry = [.2 .7]';
        % vrys=[.01 1 10]';
        % [x0, l_moments, min_res] = moments_method_2stage(data);
        % m1 = x0(1)*vry; 
        % s1 = x0(2)*vrys;
        % m2 = x0(3)*vry;
        % s2 = x0(4)*vrys;
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        l = min(data)*vry;
        N = length(vry);

        % prepare parameter seeds
        
        %get all pairs of the form [m(i),s(j)]
        %these pairs represent all possible unique 
        %parameter choices for each stochastc part of the cell
        %cycle.  
        pcomb = allcomb(m,s);
        %place paramter pairs into a cell.  The parameters choices for each part
        %are now indexed
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %get all triples of indices for the parameter 
        %choices for each part of the cycle to get all 
        %parameter choices for the entire cycle
        id = allcomb(1:length(pcomb),1:length(pcomb));
        %sort the pairs in ascending order.  
        %This equates choices of the form [i,j] and [j,i].
        id = sort(id,2);
        %remove repeats
        id = unique(id,'rows');
        %create a matrix of unique parameter choices for the cell cycle
        P = zeros(length(id)*N,5);
        for ii = 1:length(id)
            for jj=1:N
            P((ii-1)*N+jj,:) = [pcell{id(ii,1)},pcell{id(ii,2)},l(jj)];
            end
        end

        % optimize parameters
        pd=zeros(length(P),5);
        ld = NaN*ones(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        % if flag(i)=1 the pdf is approximated as a one-part stocahstic process
        % with a dterministic lag.
        flag=zeros(length(P),1);
        for i=1:length(P)  
            x0 = P(i,:);
             f=@(t,m1,s1,m2,s2,l)twostagepdf_lag(t,m1,s1,m2,s2,l,.01,10^(-6));
            [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            [l,flag(i)]=twostagepdf_lag(data,p(1),p(2),p(3),p(4),p(5),.01,10^(-6));
            l=sum(log(l));
            ld(i)=l
            %check error in Riemann sum
            if flag(i)==0
                [l]=twostagepdf_lag(data,p(1),p(2),p(3),p(4), p(5),.001,10^(-6));
                l=sum(log(l));
                ld_check(i)=l
            end
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_TWOSTAGELAG
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit three-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if threestagefit == 1
    % BEGIN FUNCTION FIT_THREESTAGE
    
        % prepare statistical parameters
        vry = [0.1 0.7]';
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        N = length(vry);
        
        % prepare parameter seeds
        pcomb = allcomb(m,s);
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        id = allcomb(1:length(pcomb),1:length(pcomb),1:length(pcomb));
        id = sort(id,2);
        id = unique(id,'rows');
        P = zeros(length(id),6);
        for ii = 1:length(id)
            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)},pcell{id(ii,3)}];
        end

        % optimize parameters
        pd=zeros(length(P),6);
        ld = NaN*ones(length(P),1);
        flag=zeros(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        for i=1:length(P)
            x0=P(i,:);
            g=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(x,m1,s1,m2,s2,m3,s3,.01);
            [p,conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            [l,hp(i),flag(i),E(i)]=convolv_3invG_nov(data,p(1),p(2),p(3),p(4),p(5),p(6),.01);
            l=sum(log(l));
            ld(i)=l
        end

        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_3invG_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),pd(i,6),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_THREESTAGE
end


end
