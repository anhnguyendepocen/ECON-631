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
prod_dedupe = unique(product_id);
prod_dedupe_mat = repmat(prod_dedupe',rows(product_id),1);
product_id_mat = repmat(product_id',rows(prod_dedupe),1)';

prod_fe = prod_dedupe_mat == product_id_mat;

BLP_rhs = horzcat(price,sugar,mushy,prod_fe);

%Now dedupe data

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% OLS
%%%%%%%%%%%%%%%%%%%%%%%%

beta_ols = inv(BLP_rhs' * BLP_rhs) * (BLP_rhs' * BLP_lhs);
resid = BLP_lhs - BLP_rhs * beta_ols;
beta_ols_se =  diag(sqrt( (1/rows(BLP_rhs)) * inv(BLP_rhs' * BLP_rhs) * (BLP_rhs' * resid * resid' * BLP_rhs) * inv(BLP_rhs' * BLP_rhs)));


%%
%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS: 
%%%%%%%%%%%%%%%%%%%%%%%%


