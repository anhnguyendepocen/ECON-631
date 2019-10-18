%%
%Created by RM on 2019.10.03 for  ECON 631 PS 1


%%%%%%%%%%%%%%%%%%%%%%%%
%   Import Data
%%%%%%%%%%%%%%%%%%%%%%%%

data = readtable('/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/cereal_data.xlsx');
product_id = table2array(data(:,1));
firm_id = table2array(data(:,2));
city = table2array(data(:,3));
year = table2array(data(:,4));
quarter = table2array(data(:,5));
market_share = table2array(data(:,6));
price = table2array(data(:,7));
sugar =  table2array(data(:,8));
mushy =  table2array(data(:,9));

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Create Instruments (Pull from Other Code)
%%%%%%%%%%%%%%%%%%%%%%%%
firm_city_date_group = horzcat(firm_id,city,year,quarter,findgroups(firm_id,city,year,quarter));

%first find sums at product and firm level; then at firm level
total_sugar_firm_city_date = accumarray(firm_city_date_group(:,5),sugar);
total_mushy_firm_city_date = accumarray(firm_city_date_group(:,5),mushy);
total_counts_firm_city_date = accumarray(firm_city_date_group(:,5),ones(rows(firm_city_date_group),1));

total_sugar_firm_city_date_expand = total_sugar_firm_city_date(firm_city_date_group(:,5),:);
total_mushy_firm_city_date_expand = total_mushy_firm_city_date(firm_city_date_group(:,5),:);
total_counts_firm_city_date_expand = total_counts_firm_city_date(firm_city_date_group(:,5),:);

denom = total_counts_firm_city_date_expand - 1;

avg_sugar_same_firm = (total_sugar_firm_city_date_expand - sugar) ...
                    ./ denom;
avg_mushy_same_firm = (total_mushy_firm_city_date_expand - mushy) ...
                    ./ denom;
 
instruments_with_six = horzcat(avg_sugar_same_firm,avg_mushy_same_firm);                
               
%remove firm six as no other products
not_firm_six = firm_id ~= 6; 

avg_sugar_same_firm_nosix = avg_sugar_same_firm(not_firm_six,:);
avg_mushy_same_firm_nosix = avg_mushy_same_firm(not_firm_six,:);
price_nosix = price(not_firm_six,:);
instruments = instruments_with_six(not_firm_six,:);

%%
%Create LHS for BLP

market_id = findgroups(city,year,quarter);
market_id_no_six = market_id(not_firm_six,:);

subs = findgroups(city,year,quarter);

sum_market_share = accumarray(subs,market_share);
subs_sum_market_share = horzcat(sum_market_share,unique(subs));
expand_sum_market_share = subs_sum_market_share(subs(:,1));

BLP_lhs = log(market_share) - log(1 - expand_sum_market_share);
BLP_lhs_nosix = BLP_lhs(not_firm_six,:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation Draws
%%%%%%%%%%%%%%%%%%%%%%%%

norm_rnd = normrnd(0,1,[10000,2]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimate Model
%%%%%%%%%%%%%%%%%%%%%%%%

%First guess: logit mean utility
sugar_nosix = sugar(not_firm_six,:);
mushy_nosix = mushy(not_firm_six,:);
BLP_rhs_nosix = horzcat(price_nosix, sugar_nosix, mushy_nosix, ...
                        ones(rows(mushy_nosix),1));

beta_ols = inv(BLP_rhs_nosix' * BLP_rhs_nosix) * (BLP_rhs_nosix' * BLP_lhs_nosix);
mean_utility = BLP_rhs_nosix * beta_ols;

%%
%Run minimzer;

theta0 = [.4 .5];
x = horzcat(sugar_nosix,mushy_nosix);
market_share_nosix = market_share(not_firm_six,:);
%tol = 10 ^ -11;
tol = 10 ^ -8;
sims = norm_rnd;
options = optimset('Display','iter');
[estimateblp] = fminsearch(@(theta)blp_gmm([theta],mean_utility,market_share_nosix,sims,x,price_nosix,instruments,tol,market_id_no_six),theta0,options);


%%
%Now Run Code To Get Other Parameters Given these Vals of sigma

    delta_curr = mean_utility;
    distance = 1;
    sigma = exp(estimateblp);
    shares = market_share_nosix;
    
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
   end;

   for_z = horzcat(x,instruments,ones(rows(instruments),1)); 
   P_Z_samefirm = for_z * inv(for_z' * for_z) ...
                    * for_z';
    X_nosix = horzcat(price_nosix,x,ones(rows(x),1));

    beta_blp = inv(X_nosix' *  P_Z_samefirm * X_nosix) ...
                            * (X_nosix' *  P_Z_samefirm * delta_curr);  
     %calc. next for se
    iv_resid_blp = delta_curr -  P_Z_samefirm * X_nosix * beta_blp;

   
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute Standard Errors
%%%%%%%%%%%%%%%%%%%%%%%%

estimateblp = [ 2.0582   -1.9037]

%attributes, including price
X_gmm = horzcat(price_nosix,sugar_nosix,mushy_nosix);
G0_beta = -X_gmm' * for_z / rows(for_z);

%sigmas: First
tol = 10^-6;
phi = 10^-6;
sigma_lt_one = exp(estimateblp) - phi .* [1,0];

    delta_curr_lt_one = mean_utility;
    distance = 1;
    shares = market_share_nosix;
    
   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma_lt_one' .* sims');
     share_numerator_hat = exp(delta_curr_lt_one .* ones(rows(delta_curr_lt_one),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
     market_sum_share_numerator_hat = zeros(rows(market_id), columns(share_numerator_hat));
 
     for i = 1: 10000;
        add_market_sum_share_numerator_hat = accumarray(market_id,share_numerator_hat(:,i));
        expand_add_market_sum_share_numerator_hat = add_market_sum_share_numerator_hat(market_id(:,1));
        market_sum_share_numerator_hat(:,i) = expand_add_market_sum_share_numerator_hat;
     end;
     
     share_hat = share_numerator_hat ./ (1 +  market_sum_share_numerator_hat);

     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr_lt_one + log(shares) - log(mean_share_hat);
     distance = sum(abs(delta_next - delta_curr_lt_one));
      
     delta_curr_lt_one = delta_next;

   end;

iv_resid_lt_one = delta_curr_lt_one -  P_Z_samefirm * X_nosix * beta_blp;
G0_lt_one = iv_resid_lt_one' * for_z;

sigma_gt_one = exp(estimateblp) + phi .* [1,0];

    delta_curr_gt_one = mean_utility;
    distance = 1;
    shares = market_share_nosix;
    
   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma_gt_one' .* sims');
     share_numerator_hat = exp(delta_curr_gt_one .* ones(rows(delta_curr_gt_one),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
    market_sum_share_numerator_hat = zeros(rows(market_id), columns(share_numerator_hat));
 
     for i = 1: 10000;
        add_market_sum_share_numerator_hat = accumarray(market_id,share_numerator_hat(:,i));
        expand_add_market_sum_share_numerator_hat = add_market_sum_share_numerator_hat(market_id(:,1));
        market_sum_share_numerator_hat(:,i) = expand_add_market_sum_share_numerator_hat;
     end;
     
     share_hat = share_numerator_hat ./ (1 +  market_sum_share_numerator_hat);
     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr_gt_one + log(shares) - log(mean_share_hat);
     distance = sum(abs(delta_next - delta_curr_gt_one));
      
     delta_curr_gt_one = delta_next;
     %distance

   end;

iv_resid_gt_one = delta_curr_gt_one -  P_Z_samefirm * X_nosix * beta_blp;
G0_gt_one = iv_resid_gt_one' * for_z;

G0_one = G0_gt_one - G0_lt_one;

%sigmas: Second 
phi = 10^-6;
sigma_lt_two = exp(estimateblp) - phi .* [0,1];

    delta_curr_lt_two = mean_utility;
    distance = 1;
    shares = market_share_nosix;
    
   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma_lt_two' .* sims');
     share_numerator_hat = exp(delta_curr_lt_two .* ones(rows(delta_curr_lt_two),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
     share_hat = share_numerator_hat ./ (1 +  share_numerator_hat);
     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr_lt_two + log(shares) - log(mean_share_hat);
     distance = sum(abs(delta_next - delta_curr_lt_two));
      
     delta_curr_lt_two = delta_next;

   end;

iv_resid_lt_two = delta_curr_lt_two -  P_Z_samefirm * X_nosix * beta_blp;
G0_lt_two = iv_resid_lt_two' * for_z;

sigma_gt_two = exp(estimateblp) + phi .* [0,1];

    delta_curr_gt_two = mean_utility;
    distance = 1;
    shares = market_share_nosix;
    
   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma_gt_two' .* sims');
     share_numerator_hat = exp(delta_curr_gt_two .* ones(rows(delta_curr_gt_two),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
     share_hat = share_numerator_hat ./ (1 +  share_numerator_hat);
     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr_gt_two + log(shares) - log(mean_share_hat);
     distance = sum(abs(delta_next - delta_curr_gt_two));
      
     delta_curr_gt_two = delta_next;

   end;

iv_resid_gt_two = delta_curr_gt_two -  P_Z_samefirm * X_nosix * beta_blp;
G0_gt_two = iv_resid_gt_two' * for_z;

G0_two = G0_gt_two - G0_lt_two;

G0_stack = horzcat(G0_beta',G0_one',G0_two');

Middle_Var = iv_resid_blp' * for_z * for_z' * iv_resid_blp

GMM_Total_Var = inv(G0_stack' * G0_stack) * G0_stack' * Middle_Var * G0_stack * inv(G0_stack' * G0_stack)
se_gmm = diag(sqrt(GMM_Total_Var / rows(iv_resid_blp) ) );