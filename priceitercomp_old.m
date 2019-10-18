function [prices_comp] = priceitercomp(alpha,beta,sigma,x,mc,prices_curr,norm_rnd,tol,weight_next);

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

    delta_s_delta_p_num_sim = s_curr_num_sim .* horzcat(sum_s_curr_num_sim',sum_s_curr_num_sim',sum_s_curr_num_sim')';
    delta_s_delta_p_sim = -alpha * delta_s_delta_p_num_sim ./ (s_curr_denom_sim .* s_curr_denom_sim);
    delta_s_delta_p = mean(delta_s_delta_p_sim,2);
    delta_s_delta_p_inv = delta_s_delta_p .^ (-1);

    prices_next = (mc' - delta_s_delta_p_inv .* s_curr)' ;
    distance = max(abs(prices_curr - prices_next));
    
    
    %prices_next = prices_curr;
    % have one firm update price
    iter = iter + 1;
    firm_chng = 1 + iter - 3 * floor(iter/3);
    %weight_next = .000005;
    prices_curr(firm_chng) = weight_next * prices_next(firm_chng) + (1 - weight_next) * prices_curr(firm_chng);
        
    
    if iter / 25000 == floor(iter / 25000);
        distance
        iter
    end;
    
    %prices_curr = .99995 * prices_curr + .0005 * prices_next;;

    %prices_curr
    %for i = 1:3;
    %    if prices_curr(i) < 0;
    %         prices_curr(i) = 0;
    %    end;
    %end;
    
end;

prices_comp = prices_curr;

end
