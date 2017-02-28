%% Rank-1 Simulation Study of Two-Way Functional Hazards
% This file compares tFHazards with an adhoc approach (piecewise constant hazard+tFSVD)
% on a rank-1 example.
% Input is survival data (with censoring)
% Output is two-way smooth hazard function (or hazard rates matrix)
% By Gen Li 4/15/2016

%% Generate true hazard function
n=50; % birth intervals
p=20; % life intervals
r=1; 

% true two-way hazard function
myfun=@(x)(1+exp(-(x-25).^2/125))/10; % hazard rate at birth time x
v=ones(1,p); % constant hazard function over lifetime
u=myfun(1:n);
Theta=u'*v; 



% plot hazard surface
figure();
mesh([1:p],[1:n],Theta);
xlim([0,p]);
ylim([0,n]);
xlabel('Lifetime','fontsize',25);
ylabel('Birth Time','fontsize',25);
zlabel('Hazard Rates','fontsize',25);
title('Two-Way Unit-Rank Hazard Function','fontsize',30);
zzz=get(gca,'zlim');
zlim(zzz);
set(gca,'fontsize',20);




    
%% simulate individual observations
nn=10000; % total number of observations for birth in (0,50) and life in (0,infty)
birthtime=unifrnd(0,50,1,nn);
censorind=binornd(1,0.1,1,nn); % censoring rate 0.1
hazardrate=myfun(birthtime);
lifetime=exprnd(1./hazardrate);

figure();
hist(lifetime)
xlabel('Lifetime','fontsize',20);
title('Histogram of Lifetime','fontsize',25);





%% set parameters
startbirth=10; % note: not from where we start generating data
endbirth=40;
num_b=50;
startlife=1;
endlife=20; % not set to be max(lifetime) b/c will have very sparse data for large lifetime 
num_l=80;

xgrid= startlife:((endlife-startlife)/num_l):endlife;
xgrid=xgrid(1:num_l)+(endlife-startlife)/(num_l*2); % lifetime grid
ygrid= startbirth:((endbirth-startbirth)/num_b):endbirth;
ygrid=ygrid(1:num_b)+(endbirth-startbirth)/(num_b*2); % birthtime grid
%% Convert survival data to matrix data
[D,T]=Surv2Mat(birthtime,lifetime,censorind,struct('startbirth',startbirth,...
    'endbirth',endbirth,'startlife',startlife,'endlife',endlife,'num_b',num_b,'num_l',num_l));



%% adhoc approach for estimating hazard surface 
[~,~,H_MLE] = SVDimpute(D./T);
% plot MLE of two-way piecewise constant hazard function
figure()
mesh(xgrid,ygrid,H_MLE);
xlim([startlife,endlife]);
ylim([startbirth,endbirth]);
xlabel('Lifetime','fontsize',20);
ylabel('Birth Time','fontsize',20);
zlabel('Hazard Rates','fontsize',20);
title('MLE of Two-Way Piecewise Constant Hazard Function','fontsize',25);
% zlim(zzz);
set(gca,'fontsize',20);


% tFSVD
r=1;
[u_FSVD,v_FSVD,~,~]=tFSVD(H_MLE,r);
H_adhoc=u_FSVD*v_FSVD';
% plot adhoc estimate of two-way hazard function
figure()
mesh(xgrid,ygrid,H_adhoc);
xlim([startlife,endlife]);
ylim([startbirth,endbirth]);
xlabel('Lifetime','fontsize',20);
ylabel('Birth Time','fontsize',20);
zlabel('Hazard Rates','fontsize',20);
title('Adhoc Estimation of Two-Way Hazard Function','fontsize',25);
zlim(zzz);
set(gca,'fontsize',20);

%% tFHazards
[u_ADMM,v_ADMM,H_ADMM]=tFHazards(D,T);
% plot tFHazards estimate of two-way hazard function
figure()
mesh(xgrid,ygrid,H_ADMM);
xlim([startlife,endlife]);
ylim([startbirth,endbirth]);
xlabel('Lifetime','fontsize',20);
ylabel('Birth Time','fontsize',20);
zlabel('Hazard Rates','fontsize',20);
title('tFHazard Estimation of Two-Way Hazard Function','fontsize',25);
zlim(zzz);
set(gca,'fontsize',20);
  
