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
tol = 10 ^ -8;
weight_next = .05;
pricesq1 = priceitercomp(alpha1,beta1,sigma1,x,mc,prices_curr,norm_rnd,tol,weight_next);

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
tol = 10 ^ -8;
weight_next = .05;

pricesq2 = priceitercomp(alpha2,beta2,sigma2,x,mc,prices_curr,norm_rnd,tol,weight_next);


%%
%for intuition: try with just sigma down.
%Higher sigma-more taste diff - more market power

pricesq2b_sigma = priceitercomp(alpha,beta,sigma2,x,mc,pricesq2,norm_rnd,tol,weight_next);
%so this LOWERS Price, consistent with above

pricesq2b_alpha = priceitercomp(alpha2,beta,sigma,x,mc,pricesq2,norm_rnd,tol,weight_next);
%so this INCREASES price, consistent with less price sensitive consumers

pricesq2b_beta = priceitercomp(alpha,beta2,sigma,x,mc,pricesq2,norm_rnd,tol,weight_next);
%everyone raised prices? consumers MORE price elastic (to switch to outside)
%more elastic-lower prices

prices_beta_test = zeros(50,4);
for i = 1:50;
       beta_test = beta1 - i * (1/100);
       prices_beta_test(i,1) = beta_test;
       prices_beta_test_temp = priceitercomp(alpha,beta_test,sigma,x,mc,pricesq2,norm_rnd,tol,weight_next);
       prices_beta_test(i,2:4) =  prices_beta_test_temp;
end;
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 3: Merger Price Effects
%%%%%%%%%%%%%%%%%%%%%%%%

ownership = [0 1 0; 1 0 0; 0 0 0];

prices_curr = pricesq2;
tol = 10 ^ -8;
weight_next = .05;

pricesq3 = priceitermerge(alpha2,beta2,sigma2,x,mc,prices_curr,norm_rnd,tol,weight_next,ownership);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 4: Merger Effects on Welfare
%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Calc CS: find logit inclusive values for each case
    random_coeffs = sigma * x' .* ones(rows(x'),rows(norm_rnd)) .* norm_rnd';
   
   %premerger
    pre_logit_incl_premerger_sim = exp( (beta * x - alpha * pricesq2)' .* ones(1,rows(norm_rnd))...
                    + random_coeffs );
    mean_logit_incl_premerger_sim = mean(pre_logit_incl_premerger_sim,2);
    logit_incl_premerger = log(sum(mean_logit_incl_premerger_sim));
    
    %postmerger
    pre_logit_incl_postmerger_sim = exp( (beta * x - alpha * pricesq3)' .* ones(1,rows(norm_rnd))...
                    + random_coeffs );
    mean_logit_incl_postmerger_sim = mean(pre_logit_incl_postmerger_sim,2);
    logit_incl_postmerger = log(sum(mean_logit_incl_postmerger_sim));
    
    cs_change = (1 / alpha2) * (logit_incl_postmerger - logit_incl_premerger);
    
%%
%Calc Producer Surplus
    %premerger
    s_curr_num_sim = exp( (beta * x - alpha * pricesq2)' .* ones(1,rows(norm_rnd))...
                    + random_coeffs );
    sum_s_curr_num_sim = sum(s_curr_num_sim);
    pre_s_curr_denom_sim = (1 + sum_s_curr_num_sim);
    s_curr_denom_sim = horzcat(pre_s_curr_denom_sim',pre_s_curr_denom_sim',pre_s_curr_denom_sim')';
    s_curr_sim = s_curr_num_sim ./ s_curr_denom_sim;
    s_curr = mean(s_curr_sim,2);
    share_premerger = s_curr;
    
    pre_producer_premerger = share_premerger .* (pricesq2' - mc');
    producer_premerger = sum(pre_producer_premerger);
    
    %postmerger
    s_curr_num_sim = exp( (beta * x - alpha * pricesq3)' .* ones(1,rows(norm_rnd))...
                    + random_coeffs );
    sum_s_curr_num_sim = sum(s_curr_num_sim);
    pre_s_curr_denom_sim = (1 + sum_s_curr_num_sim);
    s_curr_denom_sim = horzcat(pre_s_curr_denom_sim',pre_s_curr_denom_sim',pre_s_curr_denom_sim')';
    s_curr_sim = s_curr_num_sim ./ s_curr_denom_sim;
    s_curr = mean(s_curr_sim,2);
    share_postmerger = s_curr;
    
    pre_producer_postmerger = share_postmerger .* (pricesq3' - mc');    
    producer_postmerger = sum(pre_producer_postmerger);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%   Question 5: Merger Prices with Cost Reductions
%%%%%%%%%%%%%%%%%%%%%%%%

ownership = [0 1 0; 1 0 0; 0 0 0];

prices_curr = pricesq3;
tol = 10 ^ -8;
weight_next = .05;

mc_synergy = [.5 1 3];

pricesq5 = priceitermerge(alpha2,beta2,sigma2,x,mc_synergy,prices_curr,norm_rnd,tol,weight_next,ownership);
