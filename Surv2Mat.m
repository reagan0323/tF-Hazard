function [D,T]=Surv2Mat(birthtime,lifetime,censorind,paramstruct)
% This function converts standard (censored) survival data to two matrices
% suitable for tFHazards.m
% 
% Input: 
%   birthtime       n*1 numeric vector, starting time corresponding to n
%                   samples (e.g., onset time of a disease, arrival time of a call, etc)
% 
%   lifetime        n*1 numeric vector, duration for n samples (e.g., survival time, waiting time, etc)
% 
%   censorind       n*1 0/1 vector, 0=fully observed data, 1=right censored data
%
%   paramstruct
%         startbirth        default=min(birthtime), left end of birth time interval
%
%         endbirth          default=max(birthtime), right end of birth time interval 
%
%         startlife         default=min(lifetime), left end of lifetime interval
%
%         endlife           default=max(lifetime), right end of lifetime interval
%
%         num_b             default=round(n/100), total number of evenly spaced birth time intervals
%
%         num_l             default=round(n/100), total number of evenly spaced lifetime intervals
%
% Output: 
%   D           num_b*num_l count matrix, each entry is the total number of
%               events in that birth time/lifetime interval
%
%   T           num_b*num_l numeric matrix, each entry is the total
%               lifetime (summing up all samples) in that birth
%               time/lifetime interval
%
% Note: D and T are the required input for tFHazards.m
%
% Contact: Gen Li, PhD
%          Assistant Professor of Biostatistics, Columbia University
%          Email: gl2521@columbia.edu  
%
% CopyRight all reserved
% Last updated: 4/15/2016


n=length(birthtime);
n1=length(lifetime);
n2=length(censorind);
if n~=n1 || n~=n2
    error('Birth time, lifetime, and censoring index are not matched!');
end;


startbirth=min(birthtime);
endbirth=max(birthtime);
startlife=min(lifetime);
endlife=max(lifetime);
num_b=round(n/100);
num_l=round(n/100);

if nargin > 3 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'startbirth') ;    
    startbirth = getfield(paramstruct,'startbirth') ; 
  end ;
  if isfield(paramstruct,'endbirth') ;    
      endbirth = getfield(paramstruct,'endbirth') ; 
  end ;
  if isfield(paramstruct,'startlife') ;    
    startlife = getfield(paramstruct,'startlife') ; 
  end ;
  if isfield(paramstruct,'endlife') ;    
      endlife = getfield(paramstruct,'endlife') ; 
  end ;
  if isfield(paramstruct,'num_b') ;    
    num_b= getfield(paramstruct,'num_b') ; 
  end ; 
  if isfield(paramstruct,'num_l') ;    
    num_l= getfield(paramstruct,'num_l') ; 
  end ; 
end;

%%%%%%%%%%%%%
% optional manipulation and output
% any observation <startbirth or >endbirth or <startlife won't be
% considered
% any obervation >endlife will be marked as right censored
keep=(lifetime>startlife & birthtime>=startbirth & birthtime<endbirth);
newlifetime=lifetime(keep);
newbirthtime=birthtime(keep);
newcensorind=censorind(keep);
newcensorind(newlifetime>endlife)=1;  
disp(['Total number of qualified observations: ',num2str(length(newbirthtime))]);
disp(['Censoring rate: ',num2str(sum(newcensorind)/length(newcensorind))]);
%%%%%%%%%%%%%



% initialize D and T
D=zeros(num_b,num_l); 
T=D;

xgrid=startlife:((endlife-startlife)/num_l):endlife; % length num_l +1
ygrid=startbirth:((endbirth-startbirth)/num_b):endbirth; % length num_b +1


for i=1:num_b % each birthtime interval
    % find samples in that interval
    keep=(birthtime>=ygrid(i) & birthtime<ygrid(i+1)); % [left, right)
    lifetime_curr=lifetime(keep);
    censorind_curr=censorind(keep);
%     n_curr=length(lifetime_curr);
    
    for j=1:num_l % each lifetime interval
        sample1=max(lifetime_curr-xgrid(j),0);
        sample2=max(lifetime_curr-xgrid(j+1),0);
        T(i,j)=sum(sample1-sample2); % total time in this interval of this period
        D(i,j)=sum(censorind_curr==0 & lifetime_curr>xgrid(j) & lifetime_curr<=xgrid(j+1)); % total events in this (] interval
    end;
end;

if sum(sum(T==0))>0
    warning('T matrix has zero-valued entries. Be cautious when applying tFSVD.');
end;