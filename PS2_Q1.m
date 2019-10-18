%%
%Created by RM on 2019.10.18 for ECON 631 PS 2
%Compute pricing equilibria

%%

%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 0: Simulate Normal Dist. for Rand. Coefs
%%%%%%%%%%%%%%%%%%%%%%%%

norm_rnd = normrnd(0,1,[10000,3]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 1: Initialize Parameters
%%%%%%%%%%%%%%%%%%%%%%%%

alpha1 = 1;
beta1 = 1;
sigma1 = 1;
x = [1 2 3];
mc = [x(1,1) x(1,2) x(1,3)];

%%
%for debugging
alpha = 1;
beta = 1;
sigma = 1;
iter = 0;
prices_curr = [1 1 1];


%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 1: Compute prices
%%%%%%%%%%%%%%%%%%%%%%%%

prices_curr = [1 1 1];
tol = 10 ^ -6;
weight_next = .005;
pricesq1 = priceitercomp_old(alpha1,beta1,sigma1,x,mc,prices_curr,norm_rnd,tol,weight_next);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 2: Initialize Parameters
%%%%%%%%%%%%%%%%%%%%%%%%

alpha2 = .5;
beta2 = .5;
sigma2 = .5;
x = [1 2 3];
mc = [x(1,1) x(1,2) x(1,3)];

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 2: Compute prices
%%%%%%%%%%%%%%%%%%%%%%%%

prices_curr = pricesq1;
tol = 10 ^ -6;
weight_next = .000001;

pricesq2 = priceitercomp_old(alpha2,beta2,sigma2,x,mc,prices_curr,norm_rnd,tol,weight_next);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 3: Merger
%%%%%%%%%%%%%%%%%%%%%%%%

ownership = [1 1 0; 1 1 0; 0 0 1];