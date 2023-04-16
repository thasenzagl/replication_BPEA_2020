# delimit ;
clear all;
macro drop _all;
cap log close;
cap file close _all;
cap program drop _all;
set scrollbufsize 1000000;
set more off;


* Benchmark linear regression;
* MPM 2020-05-03;


*** SETTINGS ***;

* Files;
local dat_fred = "../data/us/output/data_export.dta"; // FRED data (for GDP series);
local dat_factor = "../data/us/factors/factors.csv"; // Factor data;
local tab_out = "output/fig_tab/benchmark_regr"; // Output table;

* Regression specification;
local horzs = "1 4";
local numlag = 4; // Number of lags in largest specification;
local reg_sample = "1975q2,2019q2"; // Sample;
local newey_lag = 4;

*** END SETTINGS ***;


*** LOAD DATA ***;

* Factor data;
import delim using `dat_factor', varnames(1) case(preserve);
rename ï date;
gen time = qofd(date(date, "MDY", 2019));
format time %tq;
drop date;

* Merge with GDP data;
merge 1:1 time using `dat_fred', keepusing(GDP) nogen;

* Clean up sample;
tsset time;
order time;
keep if tin(`reg_sample');
compress;

* Standardize factors;
egen global = std(Globalfactor);
egen financial = std(FinancialFactor);

* Annualize GDP growth;
replace GDP=4*GDP;


*** REGRESSIONS ***;

local lagm1 = `numlag'-1;

foreach h in `horzs' {;

	tssmooth ma Y=GDP, window(`=`h'-1' 1 0);

	newey F`h'.Y GDP global, lag(`newey_lag'); // Time series regression;
	est store reg_`h'_fin0_lag1;
	reg F`h'.Y GDP global; // Also get R-squared values (which are not computed by "newey" command);
	est store reg2_`h'_fin0_lag1;
	
	newey F`h'.Y GDP global financial, lag(`newey_lag');
	est store reg_`h'_fin1_lag1;
	reg F`h'.Y GDP global financial;
	est store reg2_`h'_fin1_lag1;
	
	newey F`h'.Y GDP L(0/`lagm1').global, lag(`newey_lag');
	est store reg_`h'_fin0_lag`numlag';
	reg F`h'.Y GDP L(0/`lagm1').global;
	est store reg2_`h'_fin0_lag`numlag';
	
	newey F`h'.Y GDP L(0/`lagm1').global L(0/`lagm1').financial, lag(`newey_lag');
	est store reg_`h'_fin1_lag`numlag';
	reg F`h'.Y GDP L(0/`lagm1').global L(0/`lagm1').financial;
	est store reg2_`h'_fin1_lag`numlag';
	
	drop Y;
	
};


** SAVE REGRESSION TABLE ***;

outreg2 [reg_*] using `tab_out', nocons excel tex replace; // Newey-West s.e.;
outreg2 [reg2_*] using `tab_out'_R2, e(r2_a) nocons excel tex replace; // R-squared;







