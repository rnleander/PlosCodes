function IMT_analysis_April2017(model)
% This file fits EMG, one-, two-, and three-stage stochasttic models to IMT data.

%clear all


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
fprintf('FUCCI Data\n\n');
%GetProcessedDataParts
%old_enviroment = load('FUCCI_April2017.mat');
%data=imt_b;
%data=old_enviroment.G2Time_b;
%data = csvread('fucci.csv');
data = [11.9000;10.6000;11.6000;9.8000;9.3000;9.0000;11.6000;11.1000;12.4000;13.7000;12.4000;11.9000;10.3000;12.9000;14.7000;11.6000;13.4000;13.4000;11.6000;10.3000;9.3000;13.7000;9.6000;10.1000;9.8000;10.9000;16.0000;9.3000;9.6000;10.3000;11.4000;10.6000;8.5000;10.3000;11.1000;8.0000;10.6000;7.5000;12.9000;9.0000;8.5000;12.4000;11.6000;9.6000;9.6000;14.7000;9.8000;10.3000;12.1000;8.8000;10.6000;12.1000;13.4000;12.4000;8.8000;13.2000;10.1000;11.6000;11.1000;15.8000;12.1000;12.7000;12.7000;11.1000;13.2000;11.9000;12.4000;13.2000;14.0000;8.0000;8.8000;9.3000;16.5000;14.5000;10.1000;14.2000;7.8000;13.2000;8.8000;8.8000;10.1000;11.9000;12.9000;14.5000;10.9000;10.6000;14.0000;8.8000;8.8000;9.0000;10.9000;14.5000;9.6000;12.4000;11.9000;12.4000;11.1000;14.5000;10.3000;12.4000;12.7000;11.9000;10.3000;13.7000;15.5000;14.5000;11.6000;10.6000;15.5000;14.7000;8.8000;11.6000;8.3000;17.6000;12.4000;11.6000;15.0000;13.7000;12.7000;10.9000;7.2000;8.5000;8.3000;9.6000;11.4000;12.9000;11.6000;13.4000;10.1000;11.6000;8.8000;12.4000;10.3000;16.3000;10.9000;10.1000;8.8000;9.3000;15.2000;8.5000;11.1000;8.3000;11.4000;11.9000;9.3000;9.8000;16.3000;12.7000;9.0000;11.9000;9.3000;10.3000;13.4000;11.4000;12.9000;12.4000;9.6000;10.3000;13.2000;10.6000;9.8000;11.9000;14.2000;13.4000;9.3000;9.6000;12.1000;11.9000;10.1000;14.0000;12.9000;21.7000;11.6000;12.1000;10.3000;9.8000;14.2000;13.7000;7.2000;10.9000;10.1000;9.6000;13.4000;13.2000;16.3000;11.6000;14.0000;10.9000;14.2000;12.4000;12.4000;13.4000;17.6000;10.1000;10.9000;14.0000;12.9000;9.0000;13.4000;15.0000;16.0000;8.0000;9.8000;12.4000;8.5000;9.6000;12.7000;12.1000;15.0000;16.0000;10.9000;14.2000;13.7000;11.9000;16.8000;11.4000;13.4000;12.4000;22.0000;12.4000;16.8000;12.1000;10.3000;13.4000;11.6000;10.1000;14.5000;10.6000;11.9000;15.5000;9.8000;12.4000;10.1000;8.0000;9.0000;9.3000;13.2000;11.1000;12.7000;12.1000;10.1000;13.2000;14.5000;10.1000;12.7000;12.9000;11.9000;12.4000;11.1000;8.5000;14.5000;16.5000;12.4000;9.0000;11.1000;9.8000;11.1000;11.1000;8.8000;13.2000;17.6000;16.8000;10.9000;12.4000;8.5000;14.7000];
%data=G1Time_b;

%for PC9 cells
%load('PC9_April2017.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose model to fit
%emgfit=0;
%twostagefit=0;
%onestagelag=0;
%onestagefit=0;
%threestagefit=0;
%twostagelag=0;
%threestagefitnomle=1;


twostagefitnomle=0;
twostagefitnomle_str='twostagefitnomle';
onestagelagnomle=0;
onestagelagnomle_str='onestagelagnomle';
onestagefitnomle=0;
onestagefitnomle_str='onestagefitnomle';
emgfitnomle=0;
emgfitnomle_str='emgfitnomle';
threestagefitnomle=0;
threestagefitnomle_str='threestage';

all_str='all';

if strcmp(model, twostagefitnomle_str)
    twostagefitnomle=1;
elseif strcmp(model, onestagelagnomle_str)
    onestagelagnomle=1;
elseif strcmp(model, onestagefitnomle_str)
    onestagefitnomle=1;
elseif strcmp(model, emgfitnomle_str)
    emgfitnomle=1;
elseif strcmp(model, threestagefitnomle_str)
    threestagefitnomle=1;
elseif strcmp(model, all_str)
    emgfitnomle=1;
    onestagefitnomle=1;
    onestagelagnomle=1;
    twostagefitnomle=1;
    threestagefitnomle=1;
else
  fprintf('No valid model was selected, quitting.\n');
end









%get sample statistics for fitting initializing the model parameters
num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));

% %Fit the EMG model with mle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% if emgfit == 1
%     % BEGIN FUNCTION FIT_EMG
%     % Input parameters: C1, C2, data
%     % Outputs: 
%         fprintf('emgfit\n\n');
%         % prepare statistical variables
%         vry = [.25 .5 .75]';  
%         c1=C1*vry;
%         c2=C2*vry;
%         %we vary the parameters so the the Gaussian and exponential parts of
%         %the cell cycle are responsible for a fraction of the total mean and
%         %variance in the IMT.
%         lam_v=1./c1;
%         mu_v=c1;
%         sig_v=c2.^.5;
%         N = length(vry);
%         
%         % prepare parameter seeds
%         pp = cell(N^3);
%         for i = 1:N
%             for j = 1:N
%                 for k = 1:N
%                     pp{N^2*(i-1)+N*(j-1)+k} = [lam_v(i), sig_v(j), mu_v(k)];
%                 end
%             end
%         end
%         P = zeros(N^3,3);
%         for ii = 1:N^3
%             P(ii,:) = pp{ii};
%         end
%         
%         % optimize parameters
%         ep = zeros(N^3,3);
%         le = -realmax*ones(N^3,1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000);
%         for i=1:length(P)
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f %f]\n',P(i,1),P(i,2),P(i,3));
%             x0=P(i,:);
%             options = statset('MaxIter',10000, 'MaxFunEvals',10000);
%             [ep(i,:),econf]=mle(data,'pdf',@emgpdf,'start',x0, 'lowerbound',[0 0 0],'upperbound',[100 100 100],'options',options);
%             fprintf('  ep=[%f %f %f]\n',ep(i,1),ep(i,2),ep(i,3));
%             econfint(:,:,i)=econf(:);
%             l=emgpdf(data,ep(i,1),ep(i,2),ep(i,3));
%             le(i)=sum(log(l));
%             fprintf('  l=%f\n\n',le(i));
%         end
%         
%         % common to each fit, consider factoring out
%         [max_le,ind_le]=max(le);
%         ep_max=ep(ind_le,:);
%         confint_max=econfint(:,:,ind_le);
%         fprintf('max_le=%f ind_le=%f\n',max_le,ind_le);
%         fprintf('ep_max=[%f %f %f]\n',ep_max(1),ep_max(2),ep_max(3));
%     % END FUNCTION FIT_EMG
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit the EMG model with fminsearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if emgfitnomle == 1
    % BEGIN FUNCTION FIT_EMG
    % Input parameters: C1, C2, data
    % Outputs: 
        fprintf('emgfitnomle\n\n');
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
            fprintf('i=%f\n',i);
            fprintf('  P=[%f %f %f]\n',P(i,1),P(i,2),P(i,3));
            x0=P(i,:);
            options = statset('MaxIter',10000, 'MaxFunEvals',10000);
            %[ep(i,:),econf]=mle(data,'pdf',@emgpdf,'start',x0, 'lowerbound',[0 0 0],'upperbound',[100 100 100],'options',options);
            f=@(x)(sum(-log(emgpdf(data,x(1),x(2),x(3)))));
            [ep(i,:),pval]=fminsearch(f, x0, options);
            fprintf('  ep=[%f %f %f]\n',ep(i,1),ep(i,2),ep(i,3));
            %econfint(:,:,i)=econf(:);
            %econfint(:,:,i)=[0.1, 0.1, 0.1, 0.1 ; 0.1, 0.1, 0.1, 0.1];
            l=emgpdf(data,ep(i,1),ep(i,2),ep(i,3));
            le(i)=sum(log(l));
            fprintf('  l=%f\n\n',le(i));
        end
        
        % common to each fit, consider factoring out
        [max_le,ind_le]=max(le);
        ep_max=ep(ind_le,:);
        %confint_max=econfint(:,:,ind_le);
        fprintf('max_le=%f ind_le=%f\n',max_le,ind_le);
        fprintf('ep_max=[%f %f %f]\n\n',ep_max(1),ep_max(2),ep_max(3));
    % END FUNCTION FIT_EMG
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Fit the one-stage model using mle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if onestagefit == 1
%     % BEGIN FUNCTION FIT_ONESTAGE
%         fprintf('onestagefit\n\n');
%         % prepare statistical variables
%         mu=1/C1;
%         sigma=(C2/C1^3)^.5;
%         vry=.5:.5:2;
%         m=mu*vry;
%         s=sigma*vry;
%         N=length(vry);
%         
%         % maybe these should be moved down into optimize parameters
%         % section?
%         pd=zeros(N^2,2);
%         ld=-realmax*ones(N^2,1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000);
%         
%         % prepare parameter seeds
%         pp = cell(N^2);
%         for i = 1:N
%             for j = 1:N
%                 pp{(i-1)*N+j} = [m(i), s(j)];
%             end
%         end
%         P = zeros(N^2,2);
%         for ii = 1:N^2
%             P(ii,:) = pp{ii};
%         end
%         
%         % optimize parameters
%         for i=1:N^2
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f]\n',P(i,1),P(i,2));
%             x0 = P(i,:);
%             [p,conf1]=mle(data,'pdf',@onestagepdf2,'start',x0, 'upperbound', [Inf Inf],'lowerbound',[0 0],'options',options);
%             fprintf('  p=[%f %f]\n',p(1),p(2));
%             pd(i,:)=p;
%             confint(:,:,i)=conf1(:);
%             l=onestagepdf2(data,p(1),p(2));
%             ld(i)=sum(log(l));
%             fprintf('  l=%f\n\n',ld(i));
%         end
%         
%         % common to each fit, consider factoring out
%         [max_ld,ind_ld]=max(ld);
%         pd_max=pd(ind_ld,:);
%         confint_max=confint(:,:,ind_ld);
%         fprintf('max_ld=%f row_ld=%f\n',max_ld,ind_ld);
%         fprintf('pd_max=[%f %f]\n',pd_max(1),pd_max(2));
%     % END FUNCTION FIT_ONESTAGE
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Fit the one-stage model using fminsearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onestagefitnomle == 1
    % BEGIN FUNCTION FIT_ONESTAGE_NOMLE
        fprintf('onestagefitnomle\n\n');
        
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
            fprintf('i=%f\n',i);
            fprintf('  P=[%f %f]\n',P(i,1),P(i,2));
            x0 = P(i,:);
            f=@(x)(sum(-log(onestagepdf2(data,x(1),x(2)))));
            [p,pval]=fminsearch(f, x0, options);
            fprintf('  p=[%f %f]\n',p(1),p(2));
            pd(i,:)=p;
            %confint(:,:,i)=conf1(:);
            %confint(:,:,i)=[0.1, 0.1, 0.1, 0.1 ; 0.1, 0.1, 0.1, 0.1];
            l=onestagepdf2(data,p(1),p(2));
            ld(i)=sum(log(l));
            fprintf('  l=%f\n\n',ld(i));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        %confint_max=confint(:,:,ind_ld);
        fprintf('max_ld=%f row_ld=%f\n',max_ld,ind_ld);
        fprintf('pd_max=[%f %f]\n\n',pd_max(1),pd_max(2));
    % END FUNCTION FIT_ONESTAGE_NOMLE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Fit one-stage model with lag with mle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if onestagelag == 1
%     % BEGIN FUNCTION FIT_ONESTAGELAG
%         fprintf('onestagefitlag\n\n');
%         % prepare statistical parameters
%         mu = C3/(3*C2^2);
%         sigma = (C3^3/(27*C2^5))^.5;
%         lag = C1-3*C2^2/C3;
%         vryv=[0.5 1 2];
%         vrym=[.25 .5 .75];
%         m = mu*vrym;
%         s = sigma*vryv;
%         lag = lag*vrym;
%         N = length(vrym);
%         
%         % prepare parameter seeds
%         pp = cell(N^3);
%         for i = 1:N
%             for j = 1:N
%                 for k = 1:N
%                 pp{N^2*(i-1)+N*(j-1)+k} = [m(i), s(j), lag(k)];
%                 end
%             end
%         end
%         P = zeros(N^3,3);
%         for ii = 1:N^3
%             P(ii,:) = pp{ii};
%         end
%         
%         % optimize parameters
%         pd = zeros(N^3,3);
%         ld = -realmax*ones(N^3,1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000);
%         for i=1:length(P)
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f %f]\n',P(i,1),P(i,2),P(i,3));
%             x0=P(i,:);
%             [p,conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options);
%             fprintf('  p=[%f %f %f]\n',p(1),p(2),p(3));
%             pd(i,:)=p;
%             confint(:,:,i)=conf1(:);
%             l=onestagepdf_lag(data,p(1),p(2),p(3));
%             ld(i)=sum(log(l));
%             fprintf('  l=%f\n\n',ld(i));
%         end
%         
%         % common to each fit, consider factoring out
%         [max_ld,ind_ld]=max(ld);
%         pd_max=pd(ind_ld,:);
%         confint_max=confint(:,:,ind_ld);
%         fprintf('max_ld=%f row_ld=%f\n',max_ld,ind_ld);
%         fprintf('pd_max=[%f %f %f]\n',pd_max(1),pd_max(2),pd_max(3));
%     % END FUNCTION FIT_ONESTAGELAG
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit one-stage model with lag with fminsearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onestagelagnomle == 1
    % BEGIN FUNCTION FIT_ONESTAGELAG
        fprintf('onestagefitlagnomle\n\n');
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
            fprintf('i=%f\n',i);
            fprintf('  P=[%f %f %f]\n',P(i,1),P(i,2),P(i,3));
            x0=P(i,:);
            %[p,conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options);
            f=@(x)(sum(-log(onestagepdf_lag(data,x(1),x(2),x(3)))));
            [p,pval]=fminsearch(f, x0, options);
            fprintf('  p=[%f %f %f]\n',p(1),p(2),p(3));
            pd(i,:)=p;
            %confint(:,:,i)=conf1(:);
            l=onestagepdf_lag(data,p(1),p(2),p(3));
            ld(i)=sum(log(l));
            fprintf('  l=%f\n\n',ld(i));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        %confint_max=confint(:,:,ind_ld);
        fprintf('max_ld=%f row_ld=%f\n',max_ld,ind_ld);
        fprintf('pd_max=[%f %f %f]\n\n',pd_max(1),pd_max(2),pd_max(3));
    % END FUNCTION FIT_ONESTAGELAG
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %Fit two-stage model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if twostagefit == 1
%     % BEGIN FUNCTION FIT_TWOSTAGE
%         fprintf('twostagefit\n');
%         % prepare statistical parameters
%         vry = [.25 .5 .75]';
%         % vrys=[.01 1 10]';
%         % [x0, l_moments, min_res] = moments_method_2stage(data);
%         % m1 = x0(1)*vry; 
%         % s1 = x0(2)*vrys;
%         % m2 = x0(3)*vry;
%         % s2 = x0(4)*vrys;
%         c1 = C1*vry;
%         c2 = C2*vry;
%         m = 1./c1;
%         s = (c2./c1.^3).^0.5;
%         N = length(vry);
% 
%         % prepare parameter seeds
%         
%         %get all pairs of the form [m(i),s(j)]
%         %these pairs represent all possible unique 
%         %parameter choices for each part of the cell
%         %cycle.  
%         pcomb = allcomb(m,s);
%         %place paramter pairs into a cell.  The parameters choices for each part
%         %are now indexed
%         pcell = cell(length(pcomb),1);
%         for i = 1:length(pcomb)
%             pcell{i} = pcomb(i,:);
%         end
%         %get all pairs of indices for the parameter 
%         %choices for each part of the cycle to get all 
%         %parameter choices for the entire cycle
%         id = allcomb(1:length(pcomb),1:length(pcomb));
%         %sort the pairs in ascending order.  
%         %This equates choices of the form [i,j] and [j,i].
%         id = sort(id,2);
%         %remove repeats
%         id = unique(id,'rows');
%         %create a matrix of unique parameter choices for the cell cycle
%         P = zeros(length(id),4);
%         for ii = 1:length(id)
%            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}];
%         end
% 
%         % optimize parameters
%         pd=zeros(length(P),4);
%         ld = NaN*ones(length(P),1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
%         flag=zeros(length(id),1);
%         for i=1:length(id)  
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4));
%             x0 = P(i,:);
%             f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
%             %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt2(x,m1,s1,m2,s2,.01,4);
%             [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options);
%             fprintf('  p=[%f %f %f %f]\n',p(1),p(2),p(3),p(4));
%             pd(i,:)=p;
%             confint(:,:,i)=conf1(:);
%             %[l,flag(i)]=convolv_2invG_small_sigma_test_var(data,p(1),p(2),p(3),p(4),.01,4);
%             [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(data,p(1),p(2),p(3),p(4),.01);
%             l=sum(log(l));
%             ld(i)=l;    
%             fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i));
%         end
%         
%         % we previously optimized with a larger step size, recalculate with
%         % a smaller stepsize after the fact
%         ld_true=zeros(length(ld),1);
%         for i=1:length(ld)
%             [l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_adapt_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),.001);
%             ld_true(i)=sum(log(l));
%         end
% 
%         % common to each fit, consider factoring out
%         [max_ld,row_ld]=max(ld_true);
%         pd_max = pd(row_ld,:);
%         confint_max=confint(:,:,row_ld);
%         fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld);
%         fprintf('pd_max=[%f %f %f %f]\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4));
%     % END FUNCTION FIT_TWOSTAGE
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit two-stage model without mle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if twostagefitnomle == 1
    % BEGIN FUNCTION FIT_TWOSTAGE
        fprintf('twostagefitnomle\n');
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
        %pcomb = allcomb(m,s)
        
        pcomb = [0.3397, 0.2359; 0.3397, 0.1180; 0.3397, 0.0786; 0.1699, 0.2359; 0.1699, 0.1180; 0.1699, 0.0786; 0.1132, 0.2359; 0.1132, 0.1180; 0.1132, 0.0786];
        
        
        %place paramter pairs into a cell.  The parameters choices for each part
        %are now indexed
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %get all pairs of indices for the parameter 
        %choices for each part of the cycle to get all 
        %parameter choices for the entire cycle
        %id = allcomb(1:length(pcomb),1:length(pcomb))
        
        %id = [1,1;1,2;1,3;1,4;1,5;1,6;1,7;1,8;1,9;2,1;2,2;2,3;2,4;2,5;2,6;2,7;2,8;2,9;3,1;3,2;3,3;3,4;3,5;3,6;3,7;3,8;3,9;4,1;4,2;4,3;4,4;4,5;4,6;4,7;4,8;4,9;5,1;5,2;5,3;5,4;5,5;5,6;5,7;5,8;5,9;6,1;6,2;6,3;6,4;6,5;6,6;6,7;6,8;6,9;7,1;7,2;7,3;7,4;7,5;7,6;7,7;7,8;7,9;8,1;8,2;8,3;8,4;8,5;8,6;8,7;8,8;8,9;9,1;9,2;9,3;9,4;9,5;9,6;9,7;9,8;9,9];
        %sort the pairs in ascending order.  
        %This equates choices of the form [i,j] and [j,i].
        %id = sort(id,2);
        %remove repeats
        %id = unique(id,'rows')
        id = [1,1;1,2;1,3;1,4;1,5;1,6;1,7;1,8;1,9;2,2;2,3;2,4;2,5;2,6;2,7;2,8;2,9;3,3;3,4;3,5;3,6;3,7;3,8;3,9;4,4;4,5;4,6;4,7;4,8;4,9;5,5;5,6;5,7;5,8;5,9;6,6;6,7;6,8;6,9;7,7;7,8;7,9;8,8;8,9;9,9];
        %create a matrix of unique parameter choices for the cell cycle
        P = zeros(length(id),4);
        for ii = 1:length(id)
            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}];
        end

        confint=nan(2,4,length(id));
        hp=nan(1,length(id));
        E=nan(1,length(id));
        hp_true=nan(1,length(id));
        flag_true=nan(1,length(id));
        E_true=nan(1,length(id));
        
        % optimize parameters
        pd=zeros(length(P),4);
        ld = NaN*ones(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        flag=zeros(length(id),1);
        for i=1:length(id)  
            fprintf('i=%f\n',i);
            fprintf('  P=[%f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4));
            x0 = P(i,:);
            
            %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
            
            f=@(x)(sum(-log(convolv_2invG_adapt_nov(data,x(1),x(2),x(3),x(4),.01))));

            
            % p found best pdf parameters
            % conf1  95% confidence intervals for the parameters
            
            %[p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)
           
            %[p,pval]=fmincon(f,x0,[],[],[],[],[0 0 0 0],[.2 4 .2 4],[],options);
            
            [p,pval]=fminsearch(f, x0, options);
            
            fprintf('  p=[%f %f %f %f]\n',p(1),p(2),p(3),p(4));
            pd(i,:)=p;
            
            %confint(:,:,i)=conf1(:);
            
            confint(:,:,i)=[0.1, 0.1, 0.1, 0.1 ; 0.1, 0.1, 0.1, 0.1];

            
            %[l,flag(i)]=convolv_2invG_small_sigma_test_var(data,p(1),p(2),p(3),p(4),.01,4);
            [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(data,p(1),p(2),p(3),p(4),.01);
            l=sum(log(l));
            ld(i)=l;    
            fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i));
        end
        
        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        fprintf('recalculating canidate solutions with smaller stepsize\n');
        ld_true=zeros(length(ld),1);
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_adapt_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
        fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld);
        fprintf('pd_max=[%f %f %f %f]\n\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4));
    % END FUNCTION FIT_TWOSTAGE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Fit two-stage model with lag
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if twostagelag == 1
%     % BEGIN FUNCTION FIT_TWOSTAGELAG
%         fprintf('twostagefitlag\n');
%         % prepare statistical parameters
%         vry = [.2 .7]';
%         % vrys=[.01 1 10]';
%         % [x0, l_moments, min_res] = moments_method_2stage(data);
%         % m1 = x0(1)*vry; 
%         % s1 = x0(2)*vrys;
%         % m2 = x0(3)*vry;
%         % s2 = x0(4)*vrys;
%         c1 = C1*vry;
%         c2 = C2*vry;
%         m = 1./c1;
%         s = (c2./c1.^3).^0.5;
%         l = min(data)*vry;
%         N = length(vry);
% 
%         % prepare parameter seeds
%         
%         %get all pairs of the form [m(i),s(j)]
%         %these pairs represent all possible unique 
%         %parameter choices for each stochastc part of the cell
%         %cycle.  
%         pcomb = allcomb(m,s);
%         %place paramter pairs into a cell.  The parameters choices for each part
%         %are now indexed
%         pcell = cell(length(pcomb),1);
%         for i = 1:length(pcomb)
%             pcell{i} = pcomb(i,:);
%         end
%         %get all triples of indices for the parameter 
%         %choices for each part of the cycle to get all 
%         %parameter choices for the entire cycle
%         id = allcomb(1:length(pcomb),1:length(pcomb));
%         %sort the pairs in ascending order.  
%         %This equates choices of the form [i,j] and [j,i].
%         id = sort(id,2);
%         %remove repeats
%         id = unique(id,'rows');
%         %create a matrix of unique parameter choices for the cell cycle
%         P = zeros(length(id)*N,5);
%         for ii = 1:length(id)
%             for jj=1:N
%             P((ii-1)*N+jj,:) = [pcell{id(ii,1)},pcell{id(ii,2)},l(jj)];
%             end
%         end
%         % optimize parameters
%         pd=zeros(length(P),5);
%         ld = NaN*ones(length(P),1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000);
%         % if flag(i)=1 the pdf is approximated as a one-part stocahstic process
%         % with a dterministic lag.
%         flag=zeros(length(P),1);
%         for i=1:length(P)  
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4),P(i,5));
%             x0 = P(i,:);
%             f=@(t,m1,s1,m2,s2,l)twostagepdf_lag(t,m1,s1,m2,s2,l,.01,10^(-6));
%             [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0],'options',options)
%             fprintf('  p=[%f %f %f %f %f]\n',p(1),p(2),p(3),p(4),p(5));
%             pd(i,:)=p;
%             confint(:,:,i)=conf1(:);
%             [l,flag(i)]=twostagepdf_lag(data,p(1),p(2),p(3),p(4),p(5),.01,10^(-6));
%             l=sum(log(l));
%             ld(i)=l
%             
%             
%             %check error in Riemann sum
%             if flag(i)==0
%                 [l]=twostagepdf_lag(data,p(1),p(2),p(3),p(4), p(5),.001,10^(-6));
%                 l=sum(log(l));
%                 ld_check(i)=l
%             end
%             fprintf('  l=%f ld_check=%f flag=%f\n\n',l,ld_check(i),flag(i));
%         end
% 
%         % common to each fit, consider factoring out
%         [max_ld,row_ld]=max(ld);
%         pd_max = pd(row_ld,:);
%         confint_max=confint(:,:,row_ld);
%         fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld);
%         fprintf('pd_max=[%f %f %f %f]\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4),pd_max(5));
%     % END FUNCTION FIT_TWOSTAGELAG
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Fit three-stage model with mle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if threestagefit == 1
%     % BEGIN FUNCTION FIT_THREESTAGE
%         fprintf('threestagefit\n');
%         % prepare statistical parameters
%         vry = [0.1 0.7]';
%         c1 = C1*vry;
%         c2 = C2*vry;
%         m = 1./c1;
%         s = (c2./c1.^3).^0.5;
%         N = length(vry);
%         
%         % prepare parameter seeds
%         pcomb = allcomb(m,s);
%         pcell = cell(length(pcomb),1);
%         for i = 1:length(pcomb)
%             pcell{i} = pcomb(i,:);
%         end
%         id = allcomb(1:length(pcomb),1:length(pcomb),1:length(pcomb));
%         id = sort(id,2);
%         id = unique(id,'rows');
%         P = zeros(length(id),6);
%         for ii = 1:length(id)
%             P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)},pcell{id(ii,3)}];
%         end
% 
%         % optimize parameters
%         pd=zeros(length(P),6);
%         ld = NaN*ones(length(P),1);
%         flag=zeros(length(P),1);
%         options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
%         for i=1:length(P)
%             fprintf('i=%f\n',i);
%             fprintf('  P=[%f %f %f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6));
%             x0=P(i,:);
%             g=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(x,m1,s1,m2,s2,m3,s3,.01);
%             [p,conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options);
%             fprintf('  p=[%f %f %f %f %f %f]\n',p(1),p(2),p(3),p(4),p(5),p(6));
%             pd(i,:)=p;
%             confint(:,:,i)=conf1(:);
%             [l,hp(i),flag(i),E(i)]=convolv_3invG_nov(data,p(1),p(2),p(3),p(4),p(5),p(6),.01);
%             l=sum(log(l));
%             ld(i)=l;
%             fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i));
%         end
% 
%         % we previously optimized with a larger step size, recalculate with
%         % a smaller stepsize after the fact
%         for i=1:length(ld)
%             [l,hp_true(i),flag_true(i),E_true(i)]=convolv_3invG_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),pd(i,6),.001);
%             ld_true(i)=sum(log(l));
%         end
% 
%         % common to each fit, consider factoring out
%         [max_ld,row_ld]=max(ld_true);
%         pd_max = pd(row_ld,:);
%         confint_max=confint(:,:,row_ld);
%         fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld);
%         fprintf('pd_max=[%f %f %f %f %f %f]\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4),pd_max(5),pd_max(6));
%     % END FUNCTION FIT_THREESTAGE
% end
% 
% 
%Fit three-stage model with fminsearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if threestagefitnomle == 1
    % BEGIN FUNCTION FIT_THREESTAGE
        fprintf('threestagefitnomle\n');
        % prepare statistical parameters
        vry = [0.1 0.7]';
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        N = length(vry);
        
        % prepare parameter seeds
        %pcomb = allcomb(m,s)
        pcomb = [0.8494,0.5898;0.8494,0.0843;0.1213,0.5898;0.1213,0.0843];
        
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %id = allcomb(1:length(pcomb),1:length(pcomb),1:length(pcomb));
        %id = sort(id,2);
        %id = unique(id,'rows')
        id = [1,1,1;1,1,2;1,1,3;1,1,4;1,2,2;1,2,3;1,2,4;1,3,3;1,3,4;1,4,4;2,2,2;2,2,3;2,2,4;2,3,3;2,3,4;2,4,4;3,3,3;3,3,4;3,4,4;4,4,4];
        P = zeros(length(id),6);
        for ii = 1:length(id)
            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)},pcell{id(ii,3)}];
        end

        hp=nan(1,length(id));
        E=nan(1,length(id));
        hp_true=nan(1,length(id));
        flag_true=nan(1,length(id));
        E_true=nan(1,length(id));
        
        % optimize parameters
        pd=zeros(length(P),6);
        ld = NaN*ones(length(P),1);
        flag=zeros(length(P),1);
        %options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs','Display','iter');
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        % WARNING WARNING WARNING
        % I've skipped over P(18),P(19) and P(20) for the moment b/c P(18)
        % is pathological. this needs fixed!!!!
        % WARNING WARNING WARNING
        for i=1:length(P)-3
            fprintf('i=%f\n',i);
            fprintf('  P=[%f %f %f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6));
            x0=P(i,:);
            %g=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(x,m1,s1,m2,s2,m3,s3,.01);
            %[p,conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options);
            
            f=@(x)(sum(-log(convolv_3invG_nov(data,x(1),x(2),x(3),x(4),x(5),x(6),0.1))));
            [p,pval]=fminsearch(f, x0, options);
            
            fprintf('  p=[%f %f %f %f %f %f]\n',p(1),p(2),p(3),p(4),p(5),p(6));
            pd(i,:)=p;
            %confint(:,:,i)=conf1(:);
            [l,hp(i),flag(i),E(i)]=convolv_3invG_nov(data,p(1),p(2),p(3),p(4),p(5),p(6),.01);
            l=sum(log(l));
            ld(i)=l;
            fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i));
        end

        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        fprintf('Recalculating canidate solutions with smaller stepsize\n');
        ld_true=zeros(length(ld),1);
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_3invG_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),pd(i,6),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        %confint_max=confint(:,:,row_ld);
        fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld);
        fprintf('pd_max=[%f %f %f %f %f %f]\n\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4),pd_max(5),pd_max(6));
    % END FUNCTION FIT_THREESTAGE
end

end
