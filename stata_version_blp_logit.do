
/* STATA check for Logit Part 2 */

clear
import excel using "/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/cereal_data.xlsx", first

bys city year quarter: egen total_market_share = sum(share)

g blp_lhs = log(share) - log(1 - total_market_share)

reg blp_lhs sugar mushy price

/* Part B: Other Chars */

bys city year quarter firm_id: egen sum_sugar_same_firm = sum(sugar)
bys city year quarter firm_id: egen sum_mushy_same_firm = sum(mushy)

g ones = 1
bys city year quarter firm_id: egen counts = sum(ones)
g denom = counts - 1

g instrument_sugar_other_chars = (sum_sugar_same_firm - sugar) / denom
g instrument_mushy_other_chars = (sum_mushy_same_firm - mushy) / denom

reg price sugar mushy instrument_sugar_other_chars instrument_mushy_other_chars
predict price_hat_same_firm, xb

reg blp_lhs price_hat_same_firm sugar mushy, r 

predict resid_same_firm, resid

ivreg blp_lhs sugar mushy (price = instrument_sugar_other_chars instrument_mushy_other_chars), r


preserve
keep if firm_id != 6
keep blp_lhs mushy sugar ones price instrument_sugar_other_chars instrument_mushy_other_chars resid_same_firm
export excel using "/Users/russellmorton/Desktop/Coursework/Fall 2019/ECON 631/Problem Sets/cereal_data_SE_TEST.xlsx", first(var) replace
restore

/* Part C: Rival Chars */

bys city year quarter: egen sum_sugar_all_firm = sum(sugar)
bys city year quarter: egen sum_mushy_all_firm = sum(mushy)

bys city year quarter: egen counts_all_firm = sum(ones)
g denom_rival = counts_all_firm - counts

g instrument_sugar_rival = (sum_sugar_all_firm - sum_sugar_same_firm) / denom_rival
g instrument_mushy_rival = (sum_mushy_all_firm - sum_mushy_same_firm) / denom_rival

reg price sugar mushy instrument_sugar_rival instrument_mushy_rival
predict price_hat_rival_firm, xb

reg blp_lhs price_hat_rival_firm sugar mushy





