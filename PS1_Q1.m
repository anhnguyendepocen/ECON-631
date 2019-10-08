%%
%Created by RM on 2019.10.02 for  ECON 631 PS 1


%%%%%%%%%%%%%%%%%%%%%%%%
%   Import Data
%%%%%%%%%%%%%%%%%%%%%%%%

data = readtable('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/ps1.dat');
work = table2array(data(:,1));
age = table2array(data(:,2));
educ = table2array(data(:,3));
parenteduc = table2array(data(:,4));

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% q1.1:  Run Probit
%%%%%%%%%%%%%%%%%%%%%%%%

thetazerohat = 0;
thetaonehat = .02;
thetatwohat = 0;

thetanullhat =  [thetazerohat thetaonehat thetatwohat];

%%
options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on'); 
[estimateprobit,log_like,exitflag,output,Gradient,Hessianprobit] = fminunc(@(x)llprobit([x],work,age,educ),thetanullhat,options);

cov_Hessianprobit = inv(Hessianprobit);
std_probit = sqrt(diag(cov_Hessianprobit));

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% q1.2: Marginal Effect of Education
%%%%%%%%%%%%%%%%%%%%%%%%

educbar = mean(educ);
agebar = mean(age);
margeffect_educ = estimateprobit(1,3) * normpdf(-estimateprobit(1,1) - estimateprobit(1,2) * agebar -  estimateprobit(1,3) * educbar);


%%
%%%%%%%%%%%%%%%%%%%%%%%%
% q1.3:  Run Logit
%%%%%%%%%%%%%%%%%%%%%%%%

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on'); 
[estimatelogit,log_like,exitflag,output,Gradient,Hessianlogit] = fminunc(@(x)lllogit([x],work,age,educ),thetanullhat,options);

cov_Hessianlogit = inv(Hessianlogit);
std_logit = sqrt(diag(cov_Hessianlogit));

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% q1.6:  Run Systems Probit
%%%%%%%%%%%%%%%%%%%%%%%%
thetazerohat = 0;
thetaonehat = 0;
thetatwohat = 0;
thetathreehat = 0;
thetafourhat = 0;
thetafivehat = 0;
rhohat = .01;
sigmasquaredhat = .01;

thetanullhat =  [thetazerohat thetaonehat thetatwohat thetathreehat thetafourhat thetafivehat rhohat sigmasquaredhat];

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on'); 
[estimateprobitsys,log_like_probitsys,exitflag,output,Gradient,Hessianprobitsys] = fminunc(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

cov_Hessianprobitsys = inv(Hessianprobitsys);
std_probitsys = sqrt(diag(cov_Hessianprobitsys));


%%
%Verify not sensitive to initial guesses


thetazerohat = estimateprobit(1,1) ;
thetaonehat = estimateprobit(1,2) ;
thetatwohat = estimateprobit(1,3) ;
thetathreehat = 0;
thetafourhat = 0;
thetafivehat = 0;
rhohat = 1;
sigmasquaredhat = 2;

thetanullhat =  [thetazerohat thetaonehat thetatwohat thetathreehat thetafourhat thetafivehat rhohat sigmasquaredhat];

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on'); 
[estimateprobitsys2,log_like_probitsys2,exitflag,output,Gradient,Hessianprobitsys2] = fminunc(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

cov_Hessianprobitsys2 = inv(Hessianprobitsys2);
std_probitsys2 = sqrt(diag(cov_Hessianprobitsys2));

%%
%Verify not sensitive to initial guesses


thetazerohat = .1 ;
thetaonehat = .1 ;
thetatwohat = .1;
thetathreehat = .1;
thetafourhat = .1;
thetafivehat = .1;
rhohat = 1;
sigmasquaredhat = 2;

thetanullhat =  [thetazerohat thetaonehat thetatwohat thetathreehat thetafourhat thetafivehat rhohat sigmasquaredhat];

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on','MaxFunEvals',200000,'MaxIter',1000); 
[estimateprobitsys3,log_like_probitsys3,exitflag,output,Gradient,Hessianprobitsys2] = fminunc(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

cov_Hessianprobitsys3 = inv(Hessianprobitsys3);
std_probitsys3 = sqrt(diag(cov_Hessianprobitsys3));

%%
%now try with initial values but other optimizer
%%%%%%%%%%%%%%%%%%%%%%%%
thetazerohat = 0;
thetaonehat = 0;
thetatwohat = 0;
thetathreehat = 0;
thetafourhat = 0;
thetafivehat = 0;
rhohat = .01;
sigmasquaredhat = .01;

thetanullhat =  [thetazerohat thetaonehat thetatwohat thetathreehat thetafourhat thetafivehat rhohat sigmasquaredhat];

%options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on'); 
[estimateprobitsys_search,log_like_probitsys_search] = fminsearch(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

thetazerohat = estimateprobit(1,1) ;
thetaonehat = estimateprobit(1,2) ;
thetatwohat = estimateprobit(1,3) ;
thetathreehat = 0;
thetafourhat = 0;
thetafivehat = 0;
rhohat = 1;
sigmasquaredhat = 2;

thetanullhat =  [thetazerohat thetaonehat thetatwohat thetathreehat thetafourhat thetafivehat rhohat sigmasquaredhat];

[estimateprobitsys_search2,log_like_probitsys_search2] = fminsearch(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

thetanullhat =  estimateprobitsys2;

[estimateprobitsys_search3,log_like_probitsys_search3] = fminsearch(@(x)llprobitsys([x],work,age,educ,parenteduc),thetanullhat,options);

%%
out_probitsys2 = round(estimateprobitsys2,4);
