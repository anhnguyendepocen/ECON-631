/* 
Created by RM for PS 3 for ECON 631
2019.11.22
Q1
*/

global data "/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets"

clear
import delimited using "$data/GMdata.csv", delim(",")

*sales, employment, capital, R&Dcapital  and  investment

g sales = exp(ldsalê)
g emp = exp(lemp)
g capital = exp(ldnpt)
g rdcapital = exp(ldrst)
g inv = exp(ldinv)

g ones = 1
bys indexê: egen count_firm = sum(ones)
g bal = count_firm > 3.5

local sumvars = "sales emp capital rdcapital inv"

local i = 0

g var = ""

g var_p25 = .
g var_median = .
g var_p75 = .

g var_mean = .
g var_sd = .

g var_count = .

g var_bal_p25 = .
g var_bal_median = .
g var_bal_p75 = .

g var_bal_mean = .
g var_bal_sd = .

g var_bal_count = .
gsort -bal

g obs = [_n]

foreach var of local sumvars {
 
 local i = `i' + 1
  
 egen pre_var_p25 = pctile(`var'), p(25)
 egen pre_var_median = median(`var')
 egen pre_var_p75 = pctile(`var'), p(75)
 
 egen pre_var_mean = mean(`var')
 egen pre_var_sd = sd(`var')
 
 g `var'_ones = 1 if `var' != .
 egen pre_var_count = sum(`var'_ones)
 
 g bal_`var' = `var' if bal == 1
 
 egen pre_var_bal_p25 = pctile(bal_`var'), p(25)
 egen pre_var_bal_median = median(bal_`var')
 egen pre_var_bal_p75 = pctile(bal_`var'), p(75)
 
 egen pre_var_bal_mean = mean(bal_`var')
 egen pre_var_bal_sd = sd(bal_`var')
 
  g bal_`var'_ones = 1 if `var' != . & bal > .05
 egen pre_var_bal_count = sum( bal_`var'_ones )
 
 replace var = "`var'" if obs == `i'
 
 replace var_p25 = pre_var_p25 if obs == `i'
 replace var_median = pre_var_median if obs == `i'
 replace var_p75 = pre_var_p75 if obs == `i'

 replace var_mean = pre_var_mean if obs == `i'
 replace var_sd = pre_var_sd if obs == `i'
 
 replace var_count = pre_var_count if obs == `i'

 replace var_bal_p25 = pre_var_bal_p25 if obs == `i'
 replace var_bal_median = pre_var_bal_median if obs == `i'
 replace var_bal_p75 = pre_var_bal_p75 if obs == `i'

 replace var_bal_mean = pre_var_bal_mean if obs == `i'
 replace var_bal_sd = pre_var_bal_sd if obs == `i'
 
 replace var_bal_count = pre_var_bal_count if obs == `i'
 
 capture drop pre*
 
 }
 

keep var*
drop if var_median == .

export excel using "$data/PS3_Summary_Stats.xlsx", replace first(var)



