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
%  Simulation Draws
%%%%%%%%%%%%%%%%%%%%%%%%

norm_rnd = normrnd(0,1,[10000,2]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Prepare Instruments
%%%%%%%%%%%%%%%%%%%%%%%%

firm_city_date_group = horzcat(firm_id,city,year,quarter,findgroups(firm_id,city,year,quarter));

%first find sums at product and firm level; then at firm level
total_sugar_firm_city_date = accumarray(firm_city_date_group(:,5),sugar);
total_mushy_firm_city_date = accumarray(firm_city_date_group(:,5),mushy);
total_counts_firm_city_date = accumarray(firm_city_date_group(:,5),ones(rows(firm_city_date_group),1));

total_sugar_firm_city_date_expand = total_sugar_firm_city_date(firm_city_date_group(:,5),:);
total_mushy_firm_city_date_expand = total_mushy_firm_city_date(firm_city_date_group(:,5),:);
total_counts_firm_city_date_expand = total_counts_firm_city_date(firm_city_date_group(:,5),:);

city_date_group = horzcat(city,year,quarter,findgroups(city,year,quarter));

total_sugar_city_date = accumarray(city_date_group(:,4),sugar);
total_mushy_city_date = accumarray(city_date_group(:,4),mushy);
total_counts_city_date = accumarray(city_date_group(:,4),ones(rows(city_date_group),1));

total_sugar_city_date_expand = total_sugar_city_date(city_date_group(:,4),:);
total_mushy_city_date_expand = total_mushy_city_date(city_date_group(:,4),:);
total_counts_city_date_expand = total_counts_city_date(city_date_group(:,4),:);

rivals_sugar_city_date = total_sugar_city_date_expand - total_sugar_firm_city_date_expand;
rivals_mushy_city_date = total_mushy_city_date_expand - total_mushy_firm_city_date_expand;
rivals_denom = total_counts_city_date_expand - total_counts_firm_city_date_expand;

avg_sugar_rivals = rivals_sugar_city_date ./ rivals_denom;
avg_mushy_rivals = rivals_mushy_city_date ./ rivals_denom;

instruments = horzcat(avg_sugar_rivals,avg_mushy_rivals);


%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimate Model
%%%%%%%%%%%%%%%%%%%%%%%%

subs = findgroups(city,year,quarter);

sum_market_share = accumarray(subs,market_share);
subs_sum_market_share = horzcat(sum_market_share,unique(subs));
expand_sum_market_share = subs_sum_market_share(subs(:,1));

BLP_lhs = log(market_share) - log(1 - expand_sum_market_share);

BLP_rhs = horzcat(price,sugar,mushy,ones(rows(price),1));

%First guess: logit mean utility

beta_ols = inv(BLP_rhs' * BLP_rhs) * (BLP_rhs' * BLP_lhs);
mean_utility = BLP_rhs * beta_ols;

%%
%Run minimzer;

theta0 = [0 0];
x = horzcat(sugar,mushy);
tol = 10 ^ -11;
sims = norm_rnd;
options = optimset('Display','iter');
[estimateblp] = fminsearch(@(theta)blp_gmm([theta],mean_utility,market_share,sims,x,price,instruments,tol),theta0,options);


%%
%Now Run Code To Get Other Parameters Given these Vals of sigma

    delta_curr = mean_utility;
    distance = 1;
    iter = 0;
    sigma = estimateblp;
    shares = market_share;
    
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

   first_stage_lhs = horzcat(instruments,ones(rows(instruments),1));
   beta_price_hat = inv(first_stage_lhs' * first_stage_lhs) * (first_stage_lhs' * price);
   price_hat = first_stage_lhs * beta_price_hat;
   
   BLP_LHS = horzcat(x,price_hat,ones(rows(x),1));
   beta_blp = inv(BLP_LHS' * BLP_LHS) * (BLP_LHS' * delta_curr);


