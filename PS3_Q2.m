%%
%Created by RM on 2019.10.22 for ECON 631 PS 3
%ACF estimation
%%

data = csvread('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/GMdata.csv',1,0);

log_sales = data(:,4);
log_emp = data(:,5);
log_cap = data(:,6);
log_rdcap = data(:,7);
log_inv = data(:,9);


%%

%Step 1: Estimate measurement error
dummy = zeros(rows(data),1) + 1;

X_step_1 = horzcat(dummy,log_emp,log_cap,log_rdcap,log_inv, ...
                        log_emp .* log_emp, log_emp .* log_cap, ...
                        log_emp .* log_rdcap, log_emp .* log_inv, ...
                        log_cap .* log_cap, log_cap .* log_rdcap, ...
                        log_cap .* log_inv, log_rdcap .* log_rdcap, ...
                        log_rdcap .* log_inv, log_inv .* log_inv );
                    
beta_step_1 = inv(X_step_1' * X_step_1) * (X_step_1' * log_sales);
pred_step_1 = X_step_1 * beta_step_1;


%%

%Step 2a: Prep data for GMM estimation of structural parameters

lag_picker = [1:1:rows(log_emp)-1];
pre_lag_vars = horzcat(log_emp,log_cap,log_rdcap,pred_step_1);
lag_vars = pre_lag_vars(lag_picker,:);

subset_set = [2:1:rows(log_sales)];
subset_log_sales = log_sales(subset_set);
pre_subset_rhs = horzcat(log_emp,log_cap,log_rdcap);
subset_rhs = pre_subset_rhs(subset_set,:);

%%

%Step 2b: Run GMM estimation


%initial guesses 
x0 = [ 1 .5 .5 .5 .8];

%compute once with "bad" weights
weights = eye(4);
instruments = horzcat(zeros(rows(subset_rhs),1) + 1, ...
                       log_emp(lag_picker), log_cap(subset_set), ...
                       log_rdcap(subset_set) );

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on','MaxFunEvals',200000,'MaxIter',1000); 
[estacf] = fminunc(@(x)acfgmm(x,subset_log_sales,subset_rhs,lag_vars,instruments,weights),x0,options);

%%
%find optimal weights 
prod_acf = [estacf(1,2) estacf(1,3) estacf(1,4)];
constant = zeros(rows(subset_log_sales),1) + 1;
pseudo_resid_no_prod_est = subset_log_sales - estacf(1,1) * constant ...
                        - subset_rhs * prod_acf';
pseudo_resid_est = pseudo_resid_no_prod_est - ...
           estacf(1,5) * (lag_vars(:,4) - lag_vars(:,1:3) * prod_acf');
dim_instruments = columns(instruments);

pre_weight_matrix = zeros(dim_instruments);
for i = 1: rows(pseudo_resid_est);
    pre_weight_matrix = pre_weight_matrix + ...
      (instruments(i,:)' * pseudo_resid_est(i)) * (instruments(i,:)' * pseudo_resid_est(i))';
end;

pre_weight_matrix = pre_weight_matrix * (1 / rows(pseudo_resid_est) );
weight_matrix = inv(pre_weight_matrix);

%%
%Estimate with optimal weights
[estacf_opt] = fminunc(@(x)acfgmm(x,subset_log_sales,subset_rhs,lag_vars,instruments,weight_matrix),x0,options);

%%
%%%%%
%Now Boostrap by Creating New Bootstrap Data
%and Keeping Track of
%Estimates
%%%%%

bstrp_iters = 5000;

bstrp_matrix = zeros(bstrp_iters,length(estacf));
bstrp_matrix_opt = zeros(bstrp_iters,length(estacf_opt));

bstrp_matrix(1,:) = estacf;
bstrp_matrix_opt(1,:) = estacf_opt;

%%
obs = rows(subset_log_sales);
bstrp_id_set = ceil(random('uniform', 0, 1,[obs,bstrp_iters])*obs) ;

for j = 2: bstrp_iters;
    
    if floor(j / 10) == j / 10;
        j
    end;
    
    bstrp_id = bstrp_id_set(:,j);
    bstrp_X_strp_1= X_step_1(bstrp_id,:);
    bstrp_log_sales = log_sales(bstrp_id);
    bstrp_beta_step_1 = inv(bstrp_X_strp_1' * bstrp_X_strp_1) ...
        * (bstrp_X_strp_1' * bstrp_log_sales);
    bstrp_pred_step_1 = bstrp_X_strp_1 * bstrp_beta_step_1;
    
    bstrp_subset_log_sales = subset_log_sales(bstrp_id);
    bstrp_subset_rhs = subset_rhs(bstrp_id,:);
    bstrp_lag_vars = lag_vars(bstrp_id,:);
    bstrp_instruments = instruments(bstrp_id,:);

    options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on','MaxFunEvals',200000,'MaxIter',1000); 
    [bstrp_estacf] = fminunc(@(x)acfgmm(x,bstrp_subset_log_sales,bstrp_subset_rhs,bstrp_lag_vars,bstrp_instruments,weights),x0,options);

    bstrp_matrix(j,:) = bstrp_estacf;
    
    bstrp_prod_acf = [bstrp_estacf(1,2) bstrp_estacf(1,3) bstrp_estacf(1,4)];
    bstrp_pseudo_resid_no_prod_est = bstrp_subset_log_sales - bstrp_estacf(1,1) * constant ...
                        - bstrp_subset_rhs * bstrp_prod_acf';
    bstrp_pseudo_resid_est = bstrp_pseudo_resid_no_prod_est - ...
           bstrp_estacf(1,5) * (bstrp_lag_vars(:,4) - bstrp_lag_vars(:,1:3) * bstrp_prod_acf');

    bstrp_pre_weight_matrix = zeros(dim_instruments);
    for i = 1: rows(pseudo_resid_est);
        bstrp_pre_weight_matrix = bstrp_pre_weight_matrix + ...
          (bstrp_instruments(i,:)' * bstrp_pseudo_resid_est(i)) ...
          * (bstrp_instruments(i,:)' * bstrp_pseudo_resid_est(i))';
    end;

    bstrp_pre_weight_matrix = bstrp_pre_weight_matrix * (1 / rows(bstrp_pseudo_resid_est) );
    bstrp_weight_matrix = inv(bstrp_pre_weight_matrix);
    
    [bstrp_estacf_opt] = fminunc(@(x)acfgmm(x,bstrp_subset_log_sales,bstrp_subset_rhs,bstrp_lag_vars,bstrp_instruments,bstrp_weight_matrix),x0,options);
    bstrp_matrix_opt(j,:) = bstrp_estacf_opt;

end;

%%

bstrp_mean = mean(bstrp_matrix);
bstrp_mean_opt = mean(bstrp_matrix_opt);

bstrp_demeaned_sq = (bstrp_matrix - bstrp_mean) .^ 2;
bstrp_demeaned_sq_opt = (bstrp_matrix_opt - bstrp_mean_opt) .^ 2;

bstrp_se = sqrt((sum(bstrp_demeaned_sq) * (1/(bstrp_iters - 1))));
bstrp_opt_se = sqrt((sum(bstrp_demeaned_sq_opt) * (1/(bstrp_iters - 1))));

%%
%Final Answers

acf_final = horzcat(estacf',bstrp_se');

acf_opt_final = horzcat(estacf_opt',bstrp_opt_se');

xlswrite('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/PS3_Q2_output_I_weights.xls',acf_final)
xlswrite('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/PS3_Q2_output_opt_weights.xls',acf_opt_final)
