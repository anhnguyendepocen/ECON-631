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

%create product fe
%prod_dedupe = unique(product_id);
%prod_dedupe_mat = repmat(prod_dedupe',rows(product_id),1);
%product_id_mat = repmat(product_id',rows(prod_dedupe),1)';

%prod_fe = prod_dedupe_mat == product_id_mat;

%BLP_rhs = horzcat(price,sugar,mushy,prod_fe);
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
% 2SLS: 
%%%%%%%%%%%%%%%%%%%%%%%%
%remove firm six as no other products

not_firm_six = firm_id ~= 6; 

check_subset = firmid(not_firm_six,:);

firm_prod_city_date_group = horzcat(firm_id,product_id,city,year,quarter,findgroups(firm_id,product_id,city,year,quarter));

%first find sums at product and firm level; then at firm level
total_sugar_firm_prod_city_date = accumarray(firm_prod_group(:,6),sugar);
total_mushy_firm_prod_city_date = accumarray(firm_prod_group(:,6),mushy);
total_counts_firm_prod_city_date = accumarray(firm_prod_group(:,6),ones(rows(firm_prod_city_date_group),1));

total_sugar_firm_prod_city_date_expand = total_sugar_firm_prod(firm_prod_group(:,6),:);
total_mushy_firm_prod_city_date_expand = total_mushy_firm_prod(firm_prod_group(:,6),:);
total_counts_firm_prod_city_date_expand = total_counts_firm_prod(firm_prod_group(:,6),:);

denom = total_counts_firm_expand - total_counts_firm_prod_expand;

avg_price_same_firm = (total_price_firm_expand - total_price_firm_prod_expand) ...
                   ./ (total_counts_firm_expand - total_counts_firm_prod_expand);
avg_sugar_same_firm = (total_sugar_firm_expand - total_sugar_firm_prod_expand) ...
                    ./ (total_counts_firm_expand - total_counts_firm_prod_expand);
avg_mushy_same_firm = (total_mushy_firm_expand - total_mushy_firm_prod_expand) ...
                    ./ (total_counts_firm_expand - total_counts_firm_prod_expand);
                 
                 




