function [u,v,impX] = SVDimpute(X)
% This function use rank-1 SVD power algorithm to impute NaN in X
%
% input: 
%   X       n*p matrix, with NaN values
%
% Output: 
%   impX    n*p matrix, all real values
%   u       n*1 vector
%   v       p*1 unit-norm vector
% 
% Note: u and v are NOT svd loadings of impX, but rather of X
%
% Contact: Gen Li, PhD
%          Assistant Professor of Biostatistics, Columbia University
%          Email: gl2521@columbia.edu  
%
% CopyRight all reserved
% Last updated: 4/15/2016

if sum(isnan(X(:)))==0
    impX=X;
    [iu,id,iv]=svds(impX,1);
    u=iu*id;
    v=iv;
    return;
end;

[n,p]=size(X);
index=~isnan(X); % 1= real value

% initial 
u=rand(n,1);
u_old=zeros(n,1);
v=rand(p,1);
v_old=zeros(p,1);

% stopping rule
thres=1e-3;
niter=1000;
iter=1;
% iterative algorithm for estimating u and v
while norm(u-u_old)+norm(v-v_old)>thres & iter<niter
    u_old=u;
    v_old=v;
    
    % estimate v
    xxinv=zeros(p,1);
    xy=zeros(p,1);
    for j=1:p
        tempind=index(:,j);
        xxinv(j)= 1/(u(tempind)'*u(tempind));
        xy(j)=u(tempind)'*X(tempind,j);
    end;
    v=xxinv.*xy;    
        
    % normalize v
    v=v/norm(v);
    
    % estimate u
    xxinv=zeros(n,1);
    xy=zeros(n,1);
    for i=1:n
        tempind=index(i,:);
        xxinv(i)= 1/(v(tempind)'*v(tempind));
        xy(i)=X(i,tempind)*v(tempind);
    end;
    u=xxinv.*xy;
    
    iter=iter+1;
end;

% impute NaN values in X
r1approx=u*v';
impX=X;
impX(isnan(X))=r1approx(isnan(X));
