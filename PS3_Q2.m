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
%Step 1a: Estimate measurement error
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
%Step 1b: Find which firms have lag variables and can be used 

max_rows = rows(data);
lag_picker = zeros(rows(data),1);
nonlag_picker = zeros(rows(data),1);

for i = 2:max_rows;
    i_lag = i - 1;
    if data(i_lag,1) == data(i,1);
        lag_picker(i,1) = i_lag;
        nonlag_picker(i,1) = i;
    end;
end;

lag_picker(lag_picker == 0) = [];
nonlag_picker(nonlag_picker == 0) = [];

%%
%Step 2a: Prep data for GMM estimation of structural parameters

pre_lag_vars = horzcat(log_emp,log_cap,log_rdcap,pred_step_1);
lag_vars = pre_lag_vars(lag_picker,:);

subset_log_sales = log_sales(nonlag_picker);
pre_subset_rhs = horzcat(log_emp,log_cap,log_rdcap);
subset_rhs = pre_subset_rhs(nonlag_picker,:);

%%
%Step 2b: Run GMM estimation

%initial guesses 
%beta_constant,beta_emp,beta_cap,beta_rdcap, rho 
x0 = [.5 .5 .5 .5 .5];

%compute once with "bad" weights
weights = eye(4);
instruments = horzcat(zeros(rows(subset_rhs),1) + 1, ...
                       log_emp(lag_picker), log_cap(nonlag_picker), ...
                       log_rdcap(nonlag_picker) );

options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-20,'TolX',1e-20,'Diagnostics','on','MaxFunEvals',500000,'MaxIter',5000); 
[estacf] = fminunc(@(x)acfgmm(x,subset_log_sales,subset_rhs,lag_vars,instruments,weights),x0,options);

%%
%find optimal weights 
prod_acf = [estacf(1,2) estacf(1,3) estacf(1,4)];
constant = zeros(rows(subset_log_sales),1) + 1;
pseudo_resid_no_prod_est = subset_log_sales - estacf(1,1) * constant ...
                        - subset_rhs * prod_acf';
pseudo_resid_est = pseudo_resid_no_prod_est - ...
           estacf(1,5) * (lag_vars(:,4) - estacf(1,1) * constant - ...
                         lag_vars(:,1:3) * prod_acf');
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
x1 = estacf / 2;
[estacf_opt] = fminunc(@(x)acfgmm(x,subset_log_sales,subset_rhs,lag_vars,instruments,weight_matrix),x1,options);

%%
%%%%%
%Now Boostrap by Creating New Bootstrap Data
%and Keeping Track of
%Estimates
%%%%%

bstrp_iters = 1000;

bstrp_matrix = zeros(bstrp_iters,length(estacf));
bstrp_matrix_opt = zeros(bstrp_iters,length(estacf_opt));

bstrp_matrix(1,:) = estacf;
bstrp_matrix_opt(1,:) = estacf_opt;

%%
obs = rows(subset_log_sales);
bstrp_id_set = ceil(random('uniform', 0, 1,[obs,bstrp_iters])*obs) ;

for j = 2: bstrp_iters;

    bstrp_id = bstrp_id_set(:,j);
    bstrp_X_step_1= X_step_1(bstrp_id,:);
    bstrp_log_sales = log_sales(bstrp_id);
    bstrp_beta_step_1 = inv(bstrp_X_step_1' * bstrp_X_step_1) ...
        * (bstrp_X_step_1' * bstrp_log_sales);
    bstrp_pred_step_1 = bstrp_X_step_1 * bstrp_beta_step_1;
    
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
           bstrp_estacf(1,5) * (bstrp_lag_vars(:,4) - bstrp_estacf(1,1) * constant ...
                                            - bstrp_lag_vars(:,1:3) * bstrp_prod_acf');

    bstrp_pre_weight_matrix = zeros(dim_instruments);
    for i = 1: rows(pseudo_resid_est);
        bstrp_pre_weight_matrix = bstrp_pre_weight_matrix + ...
          (bstrp_instruments(i,:)' * bstrp_pseudo_resid_est(i)) ...
          * (bstrp_instruments(i,:)' * bstrp_pseudo_resid_est(i))';
    end;

    bstrp_pre_weight_matrix = bstrp_pre_weight_matrix * (1 / rows(bstrp_pseudo_resid_est) );
    bstrp_weight_matrix = inv(bstrp_pre_weight_matrix);
    
    x1 = bstrp_estacf / 2;
    [bstrp_estacf_opt] = fminunc(@(x)acfgmm(x,bstrp_subset_log_sales,bstrp_subset_rhs,bstrp_lag_vars,bstrp_instruments,bstrp_weight_matrix),x1,options);
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
%Block Bootstrap
block_bstrp_iters = 200;

firms_in_block = unique(data(:,1));

block_bstrp_matrix = zeros(block_bstrp_iters,length(estacf));
block_bstrp_matrix_opt = zeros(block_bstrp_iters,length(estacf_opt));

block_bstrp_matrix(1,:) = estacf;
block_bstrp_matrix_opt(1,:) = estacf_opt;

block_bstrp_id_set = ceil(random('uniform', 0, 1,[rows(firms_in_block),bstrp_iters])*rows(firms_in_block)) ;

data_as_table = array2table(data);

for j = 2: block_bstrp_iters;

    block_bstrp_id = block_bstrp_id_set(:,j);    
    firms_as_table = array2table(block_bstrp_id);
    
    block_bstrp_data = table2array(outerjoin(firms_as_table,data_as_table,'LeftKeys',1,'RightKeys',1));
    keep = 1-isnan(block_bstrp_data(:,1));
    keep_counter = (1:1:rows(keep))';
    keeping_rows = keep .* keep_counter;
    keeping_rows(keeping_rows == 0) = [];
    
    block_bstrp_data = block_bstrp_data(keeping_rows,:);
    
    block_bstrp_lag_picker = zeros(rows(block_bstrp_data),1);
    block_bstrp_nonlag_picker = zeros(rows(block_bstrp_data),1);

    for i = 2:rows(block_bstrp_data);
    i_lag = i - 1;
    if block_bstrp_data(i_lag,1) == block_bstrp_data(i,1) & block_bstrp_data(i,4) > block_bstrp_data(i_lag,4) ;
        block_bstrp_lag_picker(i,1) = i_lag;
        block_bstrp_nonlag_picker(i,1) = i;
        end;
    end;

    block_bstrp_lag_picker(block_bstrp_lag_picker == 0) = [];
    block_bstrp_nonlag_picker(block_bstrp_nonlag_picker == 0) = [];
    
    block_bstrp_log_sales = block_bstrp_data(:,5);
    log_emp = block_bstrp_data(:,6);
    log_cap = block_bstrp_data(:,7);
    log_rdcap = block_bstrp_data(:,8);
    log_inv = block_bstrp_data(:,10);
    
    block_bstrp_X_step_1 = horzcat(zeros(rows(block_bstrp_data),1) + 1,...
                        log_emp,log_cap,log_rdcap,log_inv, ...
                        log_emp .* log_emp, log_emp .* log_cap, ...
                        log_emp .* log_rdcap, log_emp .* log_inv, ...
                        log_cap .* log_cap, log_cap .* log_rdcap, ...
                        log_cap .* log_inv, log_rdcap .* log_rdcap, ...
                        log_rdcap .* log_inv, log_inv .* log_inv );
    
    
    block_bstrp_beta_step_1 = inv(block_bstrp_X_step_1' * block_bstrp_X_step_1) ...
        * (block_bstrp_X_step_1' * block_bstrp_log_sales);
    block_bstrp_pred_step_1 = block_bstrp_X_step_1 * block_bstrp_beta_step_1 ;
    
    block_bstrp_subset_log_sales = block_bstrp_log_sales(block_bstrp_nonlag_picker);
    block_bstrp_pre_subset_rhs = horzcat(log_emp,log_cap,log_rdcap);
    block_bstrp_subset_rhs = block_bstrp_pre_subset_rhs(block_bstrp_nonlag_picker,:);
    
    pre_lag_vars = horzcat(log_emp,log_cap,log_rdcap,block_bstrp_pred_step_1);
    block_bstrp_lag_vars = pre_lag_vars(block_bstrp_lag_picker,:);
    block_bstrp_instruments = horzcat(zeros(rows(block_bstrp_subset_rhs),1) + 1, ...
                       log_emp(block_bstrp_lag_picker), log_cap(block_bstrp_nonlag_picker), ...
                       log_rdcap(block_bstrp_nonlag_picker) ); 
    options  =  optimset('GradObj','off','LargeScale','off','Display','iter','TolFun',1e-14,'TolX',1e-14,'Diagnostics','on','MaxFunEvals',200000,'MaxIter',1000); 
    [block_bstrp_estacf] = fminunc(@(x)acfgmm(x,block_bstrp_subset_log_sales,block_bstrp_subset_rhs,block_bstrp_lag_vars,block_bstrp_instruments,weights),x0,options);

    block_bstrp_matrix(j,:) = block_bstrp_estacf;
    
    block_bstrp_prod_acf = [block_bstrp_estacf(1,2) block_bstrp_estacf(1,3) block_bstrp_estacf(1,4)];
    block_bstrp_pseudo_resid_no_prod_est = block_bstrp_subset_log_sales - block_bstrp_estacf(1,1) * (zeros(rows(block_bstrp_subset_log_sales)) + 1) ...
                        - block_bstrp_subset_rhs * block_bstrp_prod_acf';
    block_bstrp_pseudo_resid_est = block_bstrp_pseudo_resid_no_prod_est - ...
           block_bstrp_estacf(1,5) * (block_bstrp_lag_vars(:,4) - block_bstrp_estacf(1,1) * (zeros(rows(block_bstrp_subset_log_sales)) + 1) ...
                                            - block_bstrp_lag_vars(:,1:3) * block_bstrp_prod_acf');

    block_bstrp_pre_weight_matrix = zeros(dim_instruments);
    for i = 1: rows(block_bstrp_pseudo_resid_est);
        block_bstrp_pre_weight_matrix = block_bstrp_pre_weight_matrix + ...
          (block_bstrp_instruments(i,:)' * block_bstrp_pseudo_resid_est(i)) ...
          * (block_bstrp_instruments(i,:)' * block_bstrp_pseudo_resid_est(i))';
    end;

    block_bstrp_pre_weight_matrix = block_bstrp_pre_weight_matrix * (1 / rows(block_bstrp_pseudo_resid_est) );
    block_bstrp_weight_matrix = inv(block_bstrp_pre_weight_matrix);
    
    x1 = block_bstrp_estacf / 2;
    [block_bstrp_estacf_opt] = fminunc(@(x)acfgmm(x,block_bstrp_subset_log_sales,block_bstrp_subset_rhs,block_bstrp_lag_vars,block_bstrp_instruments,block_bstrp_weight_matrix),x1,options);
    block_bstrp_matrix_opt(j,:) = block_bstrp_estacf_opt;

end;

%%

block_bstrp_mean = mean(block_bstrp_matrix);
block_bstrp_mean_opt = mean(block_bstrp_matrix_opt);

block_bstrp_demeaned_sq = (block_bstrp_matrix - block_bstrp_mean) .^ 2;
block_bstrp_demeaned_sq_opt = (block_bstrp_matrix_opt - block_bstrp_mean_opt) .^ 2;

block_bstrp_se = sqrt((sum(block_bstrp_demeaned_sq) * (1/(block_bstrp_iters - 1))));
block_bstrp_opt_se = sqrt((sum(block_bstrp_demeaned_sq_opt) * (1/(block_bstrp_iters - 1))));



%%
%Final Answers

acf_final = horzcat(estacf',bstrp_se',block_bstrp_se');
acf_final = round(acf_final,4);

acf_opt_final = horzcat(estacf_opt',bstrp_opt_se',block_bstrp_opt_se');
acf_opt_final = round(acf_opt_final,4);

xlswrite('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/PS3_Q2_output_I_weights.xls',acf_final)
xlswrite('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/PS3_Q2_output_opt_weights.xls',acf_opt_final)


