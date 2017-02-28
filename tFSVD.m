function [U,V,Alphau,Alphav]=tFSVD(X,r,paramstruct)
% This function conducts two-way functional SVD for an observed data matrix 
% X. It will return smooth low-rank decomposition of the original matrix.
% See more details in the 2009 JASA T&M paper 
% "The Analysis of Two-Way Functional Data Using Two-Way Regularized
% Singular Value Decompositions" by Huang, Shen, Buja.
% 
% Note: sequential method, rank-specific smoothing parameters, GCV for
% tuning selection, no strict orthogonality, even-spaced discretization
%
% Input:
%   X       n*p matrix, since we introduce smoothness in both directions,
%           we don't need to specify which is sample which is variable.
%           Traditionally, rows are samples and columns are variables.
%
%   r       positive scalar, prespecified rank, the number of components to
%           be estimated using the method. Must be smaller than n and p
% 
%   paramstruct
%       Tol     scalar (default=1E-3), stopping rule for the alternating
%               algorithm for each layer
%
%       Niter   scalar (default=1E3), maximum number of iterations for the
%               alternating algorithm for each layer
%
%       Alphau  1*r numeric vector, optimal tuning parameters for U
%               If not specified (recommended), the algorithm will select 
%               the optimal tuning parameters based on GCV.
%               
%       Alphav  1*r numeric vector, optimal tuning parameters for V
%               If not specified (recommended), the algorithm will select 
%               the optimal tuning parameters based on GCV.   
%
%       AlphavRange  numeric vector (default=10.^[-3:0.01:3]), tuning range for V
%
%       AlphauRange  numeric vector (default=10.^[-3:0.01:3]), tuning range for U
%
% Output
%   U       n*r matrix, score vectors with smooth columns (no orthogonality)
%
%   V       p*r matrix, loading vectors with unit-norm smooth columns (no orthogonality)
%
%   Alphau  1*r vector, optimal tuning parameters for U. If user specify
%           Alphau in the input, this will be exactly the same as input;
%           otherwise, this is the optimal parameters selected from GCV
%
%   Alphav  1*r vector, optimal tuning parameters for V. If user specify
%           Alphau in the input, this will be exactly the same as input;
%           otherwise, this is the optimal parameters selected from GCV
%
% Contact: Gen Li, PhD
%          Assistant Professor of Biostatistics, Columbia University
%          Email: gl2521@columbia.edu  
%
% CopyRight all reserved
% Last updated: 4/15/2016

[n,p]=size(X);
if r>n || r>p
    error('Rank exceeds matrix dimension!')
end;

% default parameters
Tol=1E-3; % stopping rule
Niter=1E3; % max iterations for each layer
auflag=0; % 0: opt tuning not given
avflag=0;
Alphau=zeros(1,r); % optimal tuning for each layer by GCV
Alphav=zeros(1,r);
Alphaset=(10.^[-3:0.01:3]); % default tuning range (large enough) for all 2r parameters
Alphauset=Alphaset; % candidate for r alphau
Alphavset=Alphaset;
fhandle=0; % 0=no fig, 1=show CV plot




if nargin > 2 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'Tol') ;    
    Tol = getfield(paramstruct,'Tol') ; 
  end ;
  if isfield(paramstruct,'Niter') ;    
    Niter = getfield(paramstruct,'Niter') ; 
  end ;
  
  if isfield(paramstruct,'Fig') ;    
    fhandle= getfield(paramstruct,'Fig') ; 
  end ;
  
  if isfield(paramstruct,'Alphau') ;    %  user-specified Alphau for each layer 1*r 
    Alphau = getfield(paramstruct,'Alphau') ; 
    auflag=1; % given tunings
  end ;
  if isfield(paramstruct,'Alphav') ;    %  user-specified Alphav for each layer 1*r
    Alphav = getfield(paramstruct,'Alphav') ; 
    avflag=1;
  end ;
  if isfield(paramstruct,'AlphavRange') ;    % user-defined Alphav tuning range
    Alphavset = getfield(paramstruct,'AlphavRange') ;    
  end ;
  if isfield(paramstruct,'AlphauRange') ;    % user-defined Alphau tuning range
    Alphauset = getfield(paramstruct,'AlphauRange') ;    
  end ;
end;



% tuning range
tu=length(Alphauset);
tv=length(Alphavset);


% set Omega
Qvv=eye(p)*(-2);
Qvv=spdiags(ones(p,1),1,Qvv);
Qvv=spdiags(ones(p,1),-1,Qvv);
Qvv=Qvv(:,2:(end-1));
Rvv=eye(p-2)*(2/3);
Rvv=spdiags(ones(p-2,1)*(1/6),1,Rvv);
Rvv=spdiags(ones(p-2,1)*(1/6),-1,Rvv);
tempv=chol(inv(Rvv))*Qvv'; % use cholesky decomp to fight against matlab precision of symmetry
Omegav=tempv'*tempv; % p*p  singular matrix (rank<p)
[VOmegav,DOmegav]=eig(full(Omegav)); % convert to a full matrix before applying eig, full is faster than sparse in this case


Quu=eye(n)*(-2);
Quu=spdiags(ones(n,1),1,Quu);
Quu=spdiags(ones(n,1),-1,Quu);
Quu=Quu(:,2:(end-1));
Ruu=eye(n-2)*(2/3);
Ruu=spdiags(ones(n-2,1)*(1/6),1,Ruu);
Ruu=spdiags(ones(n-2,1)*(1/6),-1,Ruu);
tempu=chol(inv(Ruu))*Quu';
Omegau=tempu'*tempu; % n*n
[VOmegau,DOmegau]=eig(full(Omegau)); 

% sequential estimate each rank
U=zeros(n,r); % est U
V=zeros(p,r);


for k=1:r % kth layer
    % remove first k ranks
    tX=X-U*V'; % residual matrix
    
    % initial est
    [u,d,v]=svds(tX,1); % v has norm 1
    u=u*d;
    v_old=zeros(p,1);
    u_old=zeros(n,1);
    alphau=Alphau(k);
    alphav=Alphav(k);
    
    % alternating between u and v
    niter=0;
    while max(norm(u_old-u)/norm(u), norm(v_old-v)/norm(v)) > Tol && niter<Niter
        u_old=u;
        v_old=v;
        
        % estimate u and alphau
        if auflag==1 % fixed alphau
            Su=VOmegau*diag(1./(1+alphau*diag(DOmegau)))*VOmegau';
            Rv=v'*Omegav*v/(norm(v)^2);
            coeff=Su/(1+alphav*Rv);
            u=coeff*tX*v/(norm(v)^2);
        else            
            CV=zeros(1,tu);
            Urec=zeros(n,tu);
            for i=1:tu
                alphau=Alphauset(i);
                Su=VOmegau*diag(1./(1+alphau*diag(DOmegau)))*VOmegau';
                Rv=v'*Omegav*v/(norm(v)^2);
                coeff=Su/(1+alphav*Rv);
                u=coeff*tX*v/(norm(v)^2);
                CV(i)=(1/n)*(norm(tX*v/(norm(v)^2)-u,'fro')^2)/(1-(1/n)*trace(coeff))^2;
                Urec(:,i)=u;
            end;
            [~,I]=min(CV); % global CV minimum (sometimes we encouter W shape problem of tuning selection)
%             [~,tempI]=findpeaks(-[inf,CV,inf]); % local CV minimum
%             I=tempI(end)-1; % last local minimum of CV curve
            
            alphau=Alphauset(I);
            u=Urec(:,I);  
            if fhandle
              figure(123);clf;subplot(1,2,1);
              plot(log10(Alphauset),CV);
              yyy=get(gca,'ylim');
              hold on;
              plot([log10(alphau),log10(alphau)],yyy,'r-');
              xlim([min(log10(Alphauset))-0.1,max(log10(Alphauset))+0.1]);
              xlabel('log10 of Alpha U') 
            end;
                 
        end;
        
       
        

        % estimate v and alphav
        if avflag==1 % fixed alphav
            Sv=VOmegav*diag(1./(1+alphav*diag(DOmegav)))*VOmegav';
            Ru=u'*Omegau*u/(norm(u)^2);
            coeff=Sv/(1+alphau*Ru);
            v=coeff*tX'*u/(norm(u)^2);
        else 
            CV=zeros(1,tv);
            Vrec=zeros(p,tv);
            for i=1:tv % ith tuning value
                alphav=Alphavset(i);
                Sv=VOmegav*diag(1./(1+alphav*diag(DOmegav)))*VOmegav';
                Ru=u'*Omegau*u/(norm(u)^2); 
                coeff=Sv/(1+alphau*Ru);
                v=coeff*tX'*u/(norm(u)^2);
                CV(i)=(1/p)*(norm(tX'*u/norm(u)^2-v,'fro'))^2/(1-(1/p)*trace(coeff))^2;
                Vrec(:,i)=v;
            end;
            [~,I]=min(CV);   
%             [~,tempI]=findpeaks(-[inf,CV,inf]); % local CV minimum
%             I=tempI(end)-1; % last local minimum of CV curve

            alphav=Alphavset(I);
            v=Vrec(:,I);
            if fhandle
              figure(123);subplot(1,2,2);
              plot(log10(Alphavset),CV);
              yyy=get(gca,'ylim');
              hold on;
              plot([log10(alphav),log10(alphav)],yyy,'r-');
              xlim([min(log10(Alphavset))-0.1,max(log10(Alphavset))+0.1]);
              xlabel('log10 of Alpha V')    
            end;

        end;
        

        niter=niter+1;   
    end;
    
    
    % normalize u and v
    d=norm(v);
    v=v/d;
    u=u*d;
    

    if(niter==Niter)
        warning('on','all');
        warning(['FSVD (k=',num2str(k),') does NOT converge after ',...
            num2str(Niter),' iterations!']);
    else
        disp(['FSVD (k=',num2str(k),') converges after ',...
            num2str(niter),' iterations.']);      
    end;
    

    
    
    U(:,k)=u;
    V(:,k)=v;
    Alphav(k)=alphav;
    Alphau(k)=alphau;
    
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % test run this function
% n=300;p=60;
% X=sin((1:n)/3)'*cos((1:p)*pi/30)+ (1:n)'/30*((1:p)/30-1).^2+ randn(n,p);
% r=2;
% [sU,sD,sV]=svds(X,2);
% figure(1);
% subplot(1,2,1);
% plot(sU);
% subplot(1,2,2);
% plot(sV);
% 
% [U,V,Alphau,Alphav]=TwowayFSVD(X,r);
% figure(2);
% subplot(1,2,1);
% plot(U);
% subplot(1,2,2);
% plot(V);
% Alphau
% Alphav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%