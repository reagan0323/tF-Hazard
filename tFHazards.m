function [U,V,Theta]=tFHazards(D,T,paramstruct)
% This function is used to estimate a two-way low-rank smooth hazard rate 
% surface from censored data. The description of methodology can be found
% in the paper "To Wait or Not To Wait: Understanding Customer Waiting 
% Behavior in Telephone Queues" by Gen Li, Jianhua Z. Huang, and Haipeng
% Shen.
%
%
% input: 
%   D       n*p matrix, number of events in n birth time intervals and p
%           lifetime intervals.
%   T       n*p matrix, total time in each interval (time segments added from each observation)
%
%   paramstruct
%       r       positive scalar (default=1), rank of the underlying hazard matrix
%
%       rho     positive scalar (default=0.1), ADMM step size
%
%       Niter   positive scalar (default=100), max number of ADMM iterations.
%
%       FSVD_Niter      positive scalar (default=100), max number of
%                       iterations for each tFSVD procedure
%     
%       thres   positive scalar (default=0.01), converging threshold for dual and primal residuals
%
%
% Output: 
%   U       n*r matrix, orthogonal score vectors for hazard rates, 
%           corresponding to birth time. Each column has mean 1. Each entry
%           can be explained as a proportion of the baseline hazard function.
%
%   V       p*r matrix, orthogonal loading vectors of hazard rates, 
%           corresponding to lifetime. Each column can be explained as a
%           baseline hazard function.
%
%   Theta   n*p nonegative matrix, estimate of hazard rates in each
%           prespecified interval. Very close to U*V'. Smooth in both
%           directions.
%
% 
% Contact: Gen Li, PhD
%          Assistant Professor of Biostatistics, Columbia University
%          Email: gl2521@columbia.edu  
%
% CopyRight all reserved
% Last updated: 4/15/2016



% get dimension
[n,p]=size(D);

% default setting
r=1;
rho=0.1;
Niter=100;
FSVD_Niter=100; % upper limit for FSVD iterations within ADMM iterations
thres=0.01;

if nargin > 2 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'r') ;    
    r = getfield(paramstruct,'r') ; 
  end ;
  if isfield(paramstruct,'rho') ;    
    rho = getfield(paramstruct,'rho') ; 
  end ;
  if isfield(paramstruct,'Niter') ;    
    Niter = getfield(paramstruct,'Niter') ; 
  end ;
  if isfield(paramstruct,'FSVD_Niter') ;    
    FSVD_Niter = getfield(paramstruct,'FSVD_Niter') ; 
  end ;
  if isfield(paramstruct,'thres') ;    
    thres= getfield(paramstruct,'thres') ; 
  end ; 
end;



% Initial est 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,Theta] = SVDimpute(D./T);
[U,V,alphau,alphav]=tFSVD(Theta,r,struct('Niter',FSVD_Niter)); 
Lambda=zeros(n,p);

disp('------------------------------------------');
disp('ADMM initialization');
disp(['initial rho=',num2str(rho)]);

% figure(100);clf;
% subplot(1,2,1)
% mesh(1:p,1:n,Theta);
% xlabel('waiting time');
% ylabel('time-of-day');
% zlabel('Hazard MLE');
% zzz1=get(gca,'zlim');
% subplot(1,2,2)
% mesh(1:p,1:n,U*V');
% xlabel('waiting time');
% ylabel('time-of-day');
% zlabel('UV^T (FSVDH output)');
% zzz2=get(gca,'zlim');
% 
% figure(101);clf;
% subplot(1,2,1);
% plot(U);xlabel('U for Hazard');
% subplot(1,2,2);
% plot(V);xlabel('V for Hazard');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



niter=1;
diff_primal=1;
diff_dual=1;

% ADMM iterations
while max(diff_primal,diff_dual) > thres  &&  niter<Niter;
    
    disp('------------------------------------------');
    disp(['ADMM iteration: ',num2str(niter)]);

    % reset parameters
    %%%%%%%%%%%%%%%%%%%%%%%%
    U_old=U;
    V_old=V;
    Theta_old=Theta;
    Lambda_old=Lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
  
    % estimate Theta
    %%%%%%%%%%%%%%%%%%%%%%%  
    tic
    for i=1:n
        for j=1:p
            temp=T(i,j)-rho*U(i,:)*V(j,:)'+Lambda(i,j);
            Theta(i,j)=(-temp+sqrt(temp^2+4*rho*D(i,j)))/(2*rho);
        end;
    end;
    T1=toc;
    %%%%%%%%%%%%%%%%%%%%%%%%
    
  
  
    % estimate UV
    %%%%%%%%%%%%%%%%%%%%%%%%
    tic
    % FSVD with adaptive CV range     *relatively fast and stable, but may not be accurate in tuning selection*
    lowu = log10(min(alphau)); highu=log10(max(alphau));
    lowv = log10(min(alphav)); highv=log10(max(alphav));
    Alphauset=10.^(max(-4,(lowu-.5)):0.01:min((highu+.5),4));
    Alphavset=10.^(max(-4,(lowv-.5)):0.01:min((highv+.5),4));
    [U,V,alphau,alphav]=tFSVD(Theta+Lambda/rho,r,...
        struct('Niter',FSVD_Niter,'AlphauRange',Alphauset,'AlphavRange',Alphavset));    
    T2=toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    



    % update rho (a variant of traditional fixed-rho ADMM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    const1=2; % increase/decrease rho of this fold
    const2=10;
    rho_primal=norm(Theta-U*V','fro');
    rho_dual=norm(Theta-Theta_old,'fro'); % when T is large, Theta hardly changes along iterations, therefore this is small
    if rho_primal>10*rho_dual % encourage large rho to make Theta and UV similar
        rho=const2*rho; % increase rho, when primal residual is too large
    elseif rho_primal>sqrt(10)*rho_dual
        rho=const1*rho;
    elseif rho_dual>10*rho_primal % encourage small rho to stablize Theta est
        rho=rho/const2;
    elseif rho_dual>sqrt(10)*rho_primal
        rho=rho/const1;
    end;
    %check
    disp(['rho=',num2str(rho)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
       
    
    % update Lambda
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lambda=Lambda_old+rho*(Theta-U*V');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % stopping rule (relative difference of primal/dual residuals)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % residual (NOT exactly the same with Boyd's 2010)
    res_primal=norm(Theta-U*V','fro'); 
    res_dual=max(norm(U*V'-U_old*V_old','fro'),norm(Theta-Theta_old,'fro')); % from intuition
    diff_primal=res_primal/max(norm(U*V','fro'),norm(Theta,'fro'));
    diff_dual=res_dual/max(norm(U*V','fro'),norm(Theta,'fro'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % iteration output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Time(Opt, FSVD): ',num2str([T1,T2])]);
    disp(['Primal residual: ',num2str(diff_primal)]);
    disp(['Dual residual: ',num2str(diff_dual)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    niter=niter+1;
end;



% Final display
if(niter==Niter)
    warning('on','all');
    warning(['ADMM does NOT converge after ',num2str(niter),' iterations!']);
else
    disp(['ADMM converges after ',num2str(niter),' iterations.']);
end;

