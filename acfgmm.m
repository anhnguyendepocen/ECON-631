function [gmmobj] = acfgmm(x,subset_log_sales,subset_rhs,lag_vars,instruments,weights);

%x = x0;
prod_x = [x(1,2) x(1,3) x(1,4)];

constant = zeros(rows(subset_log_sales),1) + 1;
pseudo_resid_no_prod = subset_log_sales - x(1,1) * constant ...
                        - subset_rhs * prod_x';
pseudo_resid = pseudo_resid_no_prod - ...
           x(1,5) * (lag_vars(:,4) - x(1,1) * constant ... 
                                   - lag_vars(:,1:3) * prod_x');


dim_instruments = columns(instruments);
mom_distance = zeros(dim_instruments,1);

for i = 1 : rows(pseudo_resid);
    mom_distance = mom_distance + ...
       instruments(i,:)' * pseudo_resid(i);
end;

gmmobj = (mom_distance/rows(pseudo_resid) )' * weights * (mom_distance/rows(pseudo_resid) );

end
