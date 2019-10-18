function [prices_comp] = priceitermerger(alpha,beta,sigma,x,mc,ownership,prices_curr,norm_rnd,tol);

distance = 1;
iter = 0;

while distance > tol;
 
    s_curr_num_sim = exp( (beta * x - alpha * prices_curr)' .* ones(1,rows(norm_rnd))...
                    + sigma * horzcat(norm_rnd,norm_rnd,norm_rnd)' );
    sum_s_curr_num_sim = sum(s_curr_num_sim);
    pre_s_curr_denom_sim = (1 + sum_s_curr_num_sim);
    s_curr_denom_sim = horzcat(pre_s_curr_denom_sim',pre_s_curr_denom_sim',pre_s_curr_denom_sim')';
    s_curr_sim = s_curr_num_sim ./ s_curr_denom_sim;
    s_curr = mean(s_curr_sim,2);

    delta_s_delta_p_num_sim = s_curr_num_sim .* horzcat(sum_s_curr_num_sim',sum_s_curr_num_sim',sum_s_curr_num_sim')';
    delta_s_delta_p_sim = -alpha * delta_s_delta_p_num_sim ./ (s_curr_denom_sim .* s_curr_denom_sim);
    delta_s_delta_p = mean(delta_s_delta_p_sim,2);
    delta_s_delta_p_inv = delta_s_delta_p .^ (-1);
    
    delta_s_delta_p_other_sim =  alpha * sum_s_curr_num_sim ./ (s_curr_denom_sim .* s_curr_denom_sim);
    delta_s_delta_p_other = mean(delta_s_delta_p_other_sim,2);

    pre_prices_next = (mc' - delta_s_delta_p_inv .* s_curr)' ;
    prices_next = prices_curr;
    % have one firm update price
    iter = iter + 1;
    firm_chng = 1 + iter - 3 * floor(iter/3);
    prices_next(firm_chng) = .5 * pre_prices_next(firm_chng) + .5 * prices_curr(firm_chng);
    distance = max(abs(prices_curr - prices_next))
    prices_curr = prices_next;

    %for i = 1:3;
    %    if prices_curr(i) < 0;
    %         prices_curr(i) = 0;
    %    end;
    %end;
    
end;

prices_comp = prices_curr;

end
