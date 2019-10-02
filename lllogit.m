function [log_like] = lllogit(x,y,x1,x2);
    %debug 
    %x = thetanullhat;
    %y = work;
    %x1 = age;
    %x2 = educ;

    % Probability
    prob_y_1 = 1 ./ (1 + exp(-x(1,1) - x(1,2) .* x1 - x(1,3) .* x2 ) );
    
    % Log likelihood
    ln_like_vec = y .* log(prob_y_1) + (1 - y) .* log(1 - prob_y_1);
    log_like = -sum(ln_like_vec);
end