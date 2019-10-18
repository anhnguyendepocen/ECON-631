function [gmm_obj] = blp_gmm(sigma,initial,shares,sims,x,price,instruments,tol,market_id);
 
    %blp_counter = .5;
    %blp_counter = blp_counter + .5


%%
%{
    sigma = [.4 .5];
    tol = 10 ^ -8;
    shares = market_share_nosix;
    sims = norm_rnd;
    x = horzcat(sugar_nosix,mushy_nosix);
    price = price_nosix;
    initial = mean_utility;
    instruments = instruments;
    market_id = market_id_no_six;
 %}   
    %%
    %Berry Inverstion:  
    sigma_exp = exp(sigma);
    log(sigma_exp)
    
    delta_curr = initial;
    distance = 1;
    iter = 0;
    
    distance_tracker = ones(10000,1);
    
    market_id_expand = market_id .* ones(rows(delta_curr),columns(delta_curr));

   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma_exp' .* sims');
     share_numerator_hat = exp(delta_curr .* ones(rows(delta_curr),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
     market_sum_share_numerator_hat = zeros(rows(market_id), columns(share_numerator_hat));
 
     for i = 1: 10000;
        add_market_sum_share_numerator_hat = accumarray(market_id,share_numerator_hat(:,i));
        expand_add_market_sum_share_numerator_hat = add_market_sum_share_numerator_hat(market_id(:,1));
        market_sum_share_numerator_hat(:,i) = expand_add_market_sum_share_numerator_hat;
     end;
     
     share_hat = share_numerator_hat ./ (1 +  market_sum_share_numerator_hat);
     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr + log(shares) - log(mean_share_hat);
     distance = max(abs(delta_next - delta_curr));
      
     delta_curr = delta_next;
     iter = iter + 1;
     distance_tracker(iter) = distance;
     %distance
   end;
   
   %%
   %IV Regression to back out the betas and price sensitivity
   for_z = horzcat(x,instruments,ones(rows(instruments),1)); 
   P_Z_samefirm = for_z * inv(for_z' * for_z) ...
                    * for_z';
    X_nosix = horzcat(price,x,ones(rows(x),1));

    beta_2SLS_samefirm = inv(X_nosix' *  P_Z_samefirm * X_nosix) ...
                            * (X_nosix' *  P_Z_samefirm * delta_curr);  

   
   %%
   %compute obj fct val
   iv_resid =  delta_curr -  P_Z_samefirm * X_nosix * beta_2SLS_samefirm;
   pre_gmm_resid = iv_resid' * for_z;
   gmm_obj = pre_gmm_resid * pre_gmm_resid';
   
   beta_2SLS_samefirm
   gmm_obj
   %blp_counter = blp_counter + .5
   
end