function [prices_merge] = priceitermerge(alpha,beta,sigma,x,mc,prices_curr,norm_rnd,tol,weight_next,ownership);

distance = 1;
iter = 0;

while distance > tol;
 
    random_coeffs = sigma * x' .* ones(rows(x'),rows(norm_rnd)) .* norm_rnd';
    s_curr_num_sim = exp( (beta * x - alpha * prices_curr)' .* ones(1,rows(norm_rnd))...
                    + random_coeffs );
    sum_s_curr_num_sim = sum(s_curr_num_sim);
    pre_s_curr_denom_sim = (1 + sum_s_curr_num_sim);
    s_curr_denom_sim = horzcat(pre_s_curr_denom_sim',pre_s_curr_denom_sim',pre_s_curr_denom_sim')';
    s_curr_sim = s_curr_num_sim ./ s_curr_denom_sim;
    s_curr = mean(s_curr_sim,2);

    delta_s_delta_p_sim = alpha * (-s_curr_sim + ...
                            (s_curr_num_sim .* s_curr_num_sim) ./ ( s_curr_denom_sim .* s_curr_denom_sim) );
    delta_s_delta_p = mean(delta_s_delta_p_sim,2);
    delta_s_delta_p_inv = delta_s_delta_p .^ (-1);
    
    %compute derivatives to hit with ownership matrix
    delta_s_delta_p_mat = zeros(columns(x),columns(x));
    for i = 1:3;
        for j = 1:3;
            if i == j;
                delta_s_delta_p_mat(i,j) = delta_s_delta_p(i);
            else
                num_delta_s_delta_p_mat_sim = s_curr_num_sim(i,:) .*  s_curr_num_sim(j,:);
                denom_delta_s_delta_p_mat_sim = pre_s_curr_denom_sim .* pre_s_curr_denom_sim;
                delta_s_delta_p_mat(i,j) = alpha * mean( num_delta_s_delta_p_mat_sim ./ denom_delta_s_delta_p_mat_sim, 2);
            end;
        end;
    end;    
           
               
    %compute next prices
    merger_price_add = (ownership * (prices_curr - mc)') .* max(ownership .* delta_s_delta_p_mat)' ;  
    prices_next = (mc' - delta_s_delta_p_inv .* (s_curr + merger_price_add) )' ;
    distance = max(abs(prices_curr - prices_next));
    
    %prices_next = prices_curr;
    % have one firm update price
    iter = iter + 1;
    firm_chng = 1 + iter - 3 * floor(iter/3);
    prices_curr(firm_chng) = weight_next * prices_next(firm_chng) + (1 - weight_next) * prices_curr(firm_chng);
        
    
    if iter / 25000 == floor(iter / 25000);
        distance
        iter
    end;
    
end;

prices_merge = prices_curr;

end
