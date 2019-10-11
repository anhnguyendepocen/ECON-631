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
% Create Panel Data of Market Shares
%%%%%%%%%%%%%%%%%%%%%%%%

%Find Outside Option Value
subs = findgroups(city,year,quarter);

sum_market_share = accumarray(subs,market_share);
subs_sum_market_share = horzcat(sum_market_share,unique(subs));
expand_sum_market_share = subs_sum_market_share(subs(:,1));

BLP_lhs = log(market_share) - log(1 - expand_sum_market_share);

BLP_rhs = horzcat(price,sugar,mushy,ones(rows(price),1));

%Now dedupe data

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% OLS
%%%%%%%%%%%%%%%%%%%%%%%%

beta_ols = inv(BLP_rhs' * BLP_rhs) * (BLP_rhs' * BLP_lhs);
resid = BLP_lhs - BLP_rhs * beta_ols;

H0 = inv((1 / rows(BLP_rhs) ) * BLP_rhs' * BLP_rhs);

V0 = zeros(4,4);
for i = 1:rows(BLP_rhs);
    add_one = BLP_rhs(i,:)' * BLP_rhs(i,:) * resid(i)^2;
    V0 = V0 + add_one;
end;

V0 = V0 / (rows(BLP_rhs));

se_beta_ols = diag(sqrt(H0 * V0 * H0 / rows(BLP_rhs) ));


%%
%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS: OWN PRODUCTS
%%%%%%%%%%%%%%%%%%%%%%%%
%Prep Data

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
                 
%remove firm six as no other products
not_firm_six = firm_id ~= 6; 

avg_sugar_same_firm_nosix = avg_sugar_same_firm(not_firm_six,:);
avg_mushy_same_firm_nosix = avg_mushy_same_firm(not_firm_six,:);
price_nosix = price(not_firm_six,:);
BLP_lhs_nosix = BLP_lhs(not_firm_six,:);
%%
%Estimate 2SLS with Own Prods

instruments_samefirm = horzcat(avg_sugar_same_firm_nosix,avg_mushy_same_firm_nosix,...
                ones(rows(price_nosix),1));
%Fist Stage
first_stage_beta_hat_samefirm = inv(instruments_samefirm' * instruments_samefirm) * (instruments_samefirm' * price_nosix);
price_hat_samefirm = instruments_samefirm * first_stage_beta_hat_samefirm;

%Second Stage
sugar_nosix = sugar(not_firm_six,:);
mushy_nosix = mushy(not_firm_six,:);
BLP_rhs_samefirm = horzcat(price_hat_samefirm,sugar_nosix,mushy_nosix,...
                    ones(rows(sugar_nosix),1));
beta_2SLS_samefirm = inv(BLP_rhs_samefirm' * BLP_rhs_samefirm) * (BLP_rhs_samefirm' * BLP_lhs_nosix);
               
            

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS: RIVAL PRODUCTS
%%%%%%%%%%%%%%%%%%%%%%%%
%Prep Data

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

%%
%Estimate 2SLS with Rival Prods
instruments_rivalfirm = horzcat(avg_sugar_rivals,avg_mushy_rivals,...
                ones(rows(price),1));

%First Stage
first_stage_beta_hat_rival = inv(instruments_rivalfirm' * instruments_rivalfirm) * (instruments_rivalfirm' * price);
price_hat_rival = instruments_rivalfirm * first_stage_beta_hat_rival;

%Second Stage
BLP_rhs_rival = horzcat(price_hat_rival,sugar,mushy,...
                    ones(rows(sugar),1));
beta_2SLS_rival = inv(BLP_rhs_rival' * BLP_rhs_rival) * (BLP_rhs_rival' * BLP_lhs);
               
   
