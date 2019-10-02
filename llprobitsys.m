function [log_like] = llprobitsys(x,y,x1,x2,z);
    %debug 
    
    %{
    x = thetanullhat;
    y = work;
    x1 = age;
    x2 = educ;
    z = parenteduc;
    %}
    
    % Probability
    prob_y_1 = zeros(rows(y),1);
    mu = zeros(rows(y),1);
    
    for i = 1: rows(y) ;
        
    mu(i) =  (x(1,8)/x(1,7)) ...
        * (x2(i) - (x(1,4) + x(1,5) .* x1(i) + x(1,6) .* z(i) ) );
    prob_y_1(i) = 1 - normcdf(-x(1,1) - x(1,2) .* x1(i) - x(1,3) .* x2(i) , ...
        mu(i)  ...
        , 1 -  (x(1,8)^2/x(1,7) ) );
        if  prob_y_1(i) > .9999999
            prob_y_1(i) = .99;
        end;
        if  prob_y_1(i)  < ( 1 - .9999999)
            prob_y_1(i) = .01;
        end;
    end;
    
    % Log likelihood
    ln_like_vec = y .* log(prob_y_1) + (1 - y) .* log(1 - prob_y_1);
    log_like = -sum(ln_like_vec);
end