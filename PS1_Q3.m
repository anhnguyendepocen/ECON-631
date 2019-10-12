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

theta0 = [0 0];
x = horzcat(sugar_nosix,mushy_nosix);
market_share_nosix = market_share(not_firm_six,:);
tol = 10 ^ -11;
sims = norm_rnd;
options = optimset('Display','iter');
[estimateblp] = fminsearch(@(theta)blp_gmm([theta],mean_utility,market_share_nosix,sims,x,price_nosix,instruments,tol),theta0,options);


%%
%Now Run Code To Get Other Parameters Given these Vals of sigma

    delta_curr = mean_utility;
    distance = 1;
    iter = 0;
    sigma = estimateblp;
    shares = market_share_nosix;
    
    distance_tracker = ones(10000,1);

   while distance > tol;
       
     attributes = x;
     idiosyncratic_utility = attributes * (sigma' .* sims');
     share_numerator_hat = exp(delta_curr .* ones(rows(delta_curr),columns(idiosyncratic_utility)) ...
         +  idiosyncratic_utility);
     share_hat = share_numerator_hat ./ (1 +  share_numerator_hat);
     mean_share_hat = mean(share_hat,2);
     delta_next = delta_curr + log(shares) - log(mean_share_hat);
     distance = sum(abs(delta_next - delta_curr));
      
     delta_curr = delta_next;
     iter = iter + 1;
     distance_tracker(iter) = distance;
   end;

   for_z = horzcat(x,instruments); 
   P_Z_samefirm = for_z * inv(for_z' * for_z) ...
                    * for_z';
    X_nosix = horzcat(price,x,ones(rows(x),1));

    beta_blp = inv(X_nosix' *  P_Z_samefirm * X_nosix) ...
                            * (X_nosix' *  P_Z_samefirm * delta_curr);  

   

