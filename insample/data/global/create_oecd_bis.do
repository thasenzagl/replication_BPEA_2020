# delimit ;


* Create cross-country data set from OECD and BIS;
* MPM 2019-12-19;


*** SETTINGS ***;

* Files;
local file_eo = "oecd/EO106_20191218.csv"; // OECD EO;
local file_mei = "oecd/MEI_20191218.csv"; // OECD MEI;
local file_bispp = "bis/pp_selected_20191218.xlsx"; // BIS selected property prices;
local file_bistc = "bis/totcredit_20191218.xlsx"; // BIS credit;

* Currencies;
local curr = "AUD CAD CHF EUR GBP JPY SEK"; // Currency codes in EO data;

* Sample;
local end_sample = "2019q2"; // End EO sample here;

*** END SETTINGS ***;


*** OECD EO DATA ***;

* Import data;
import delim using `file_eo', varnames(1) case(preserve);
drop Frequency FREQUENCY Time Power* Reference* Flag*;
rename *LOCATION geo_str;

* Time;
gen time = quarterly(TIME, "YQ");
format time %tq;
drop TIME;

* Convert local currency variables to USD;
drop if VARIABLE == "EXCHUD";
foreach c of local curr {;
	gen EXCH = Value if VARIABLE=="EXCH" & UnitCode=="`c'";
	egen `c' = max(EXCH) if UnitCode=="`c'", by(time);
	gen sample = UnitCode=="`c'" & VARIABLE!="EXCH";
	replace Value = Value*`c' if sample==1;
	replace Unit = "US Dollar" if sample==1;
	replace UnitCode = "USD" if sample==1;
	drop EXCH `c' sample;
};

* Labels;
gen var_str = VARIABLE + "_" + UnitCode; // Variable name + unit;
gen var_lab = Variable + " " + Unit;
encode geo_str, gen(geo);
labmask geo, values(Country);
drop VARIABLE Variable Unit* Country;
rename Value value;

* Drop forecasted variables;
drop if time > quarterly("`end_sample'", "YQ");

* Sort and save;
sort geo var_str time;
order geo var_str time value;
compress;
save output/oecd_eo, replace;


*** OECD MEI DATA ***;

* Import data;
import delim using `file_mei', varnames(1) case(preserve) clear;
drop Frequency Time Unit* Power* Reference*;
rename *LOCATION geo_str;

* Time;
gen time = quarterly(TIME, "YQ") if FREQUENCY=="Q";
gen time_m = monthly(TIME, "YM") if FREQUENCY=="M";
replace time = qofd(dofm(time_m)) if FREQUENCY=="M";
format time %tq;
format time_m %tm;
drop TIME;
save output/oecd_mei, replace;

* Variable labels;
use output/oecd_eo, clear; // Import country labels from previous EO data;
collapse (first) geo_str, by(geo);
merge 1:m geo_str using output/oecd_mei, nogen; // Load MEI data again;
label val geo geo;
gen var_str = SUBJECT + "_" + MEASURE; // Variable name + unit;
gen var_lab = Subject + " " + Measure;
encode FlagCodes, gen(flag);
labmask flag, values(Flags);
encode FREQUENCY, gen(frequency);
drop SUBJECT Subject MEASURE Measure Country Flag* FREQUENCY;
rename Value value;

* Save;
sort geo var_str time time_m;
order geo var_str time time_m value;
compress;
save output/oecd_mei, replace;


*** MERGE OECD DATA SETS ***;

use output/oecd_mei, clear;

* Compute quarterly averages from monthly series;
egen count_tot = count(value), by(geo var_str time);
egen count_m = count(time_m), by(geo var_str time);
replace value = . if (count_m!=3 | count_tot==4) & frequency=="M":frequency; // Drop monthly value if not all 3 months are available, or if quarterly obs. is already available;

* Collapse to just quarterly observations, prioritizing quarterly series if available;
collapse (mean) value, by(geo var_str time geo_str var_lab);

* Append;
append using output/oecd_eo;
drop if missing(value);
replace var_lab = var_lab + " OECD"; // Include source in variable label;

* Save;
compress;
save output/global, replace;


*** BIS DATA ***;

foreach f in "`file_bispp'" "`file_bistc'" {;

	* Import data;
	import excel using `f', sheet("Quarterly Series") clear;

	* Store labels of variables;
	rename * v*;
	rename vA date;
	foreach v of varlist v* {; // Go through variables;
		replace `v' = subinstr(`v', "Q:", "", 1) if _n==4; // Drop frequency part of label;
		replace `v' = substr(`v'[_n+1], 1, 2) if _n==3; // Store country part of label in 3rd obs.;
		replace `v' = substr(`v', 4, .) if _n==4; // Drop country part of label in 4th obs.;
		replace `v' = subinstr(`v', ":", "_", .) if _n==4; // Replace : with _;
	};
	local j = 1;
	local the_vars = "";
	foreach v of varlist v* {; // Run through variables to store variable labels;
		local the_var = "`=`v'[4]'"; // Variable name;
		if strpos("`the_vars'", "`the_var'") == 0 {; // If variable has not already been stored;
			local var_`j' = "`the_var'"; // Variable name;
			local lab_`j' = "`=`v'[1]' `=`v'[2]'"; // Label;
			local j = `j' + 1;
			local the_vars = "`the_vars' `the_var'";
		};
	};
	local numlab = `j'; // Total number of variables;

	* Reshape data;
	foreach v of varlist v* {;
		rename `v' value`=`v'[3]'_`=`v'[4]';
	};
	drop if _n<=4;
	reshape long value, i(date) j(geo_var) string;
	destring, replace;
	drop if missing(value);
	gen geo_str_short = substr(geo_var, 1, 2);
	gen var_str = substr(geo_var, 4, .);

	* Time;
	gen time = qofd(date(date, "DMY"));
	format time %tq;
	drop date;

	* Country codes;
	kountry geo_str_short, from(iso2c) to(iso3c) marker;
	drop if MARKER==0;
	rename _ISO3C_ geo_str;
	drop MARKER NAMES_STD geo_var geo_str_short;
	tempfile tempf;
	save `tempf';
	use output/global, clear; // Load geographic names from previous database, keeping only OECD countries;
	collapse (first) geo_str, by(geo);
	merge 1:m geo_str using `tempf', nogen keep(match);
	label val geo geo;

	* Labels;
	gen var_lab = "";
	forvalues i=1/`numlab' {;
		replace var_lab = "`lab_`i''" if var_str == "`var_`i''";
	};
	replace var_lab = var_lab + " BIS"; // Include source in labels;

	* Append;
	append using output/global;
	save output/global, replace;

};

* Fix labels for credit variables;
gen creditpos = strpos(var_lab, "Credit");
replace var_lab = substr(var_lab, creditpos, .) if creditpos>0;
drop creditpos;


*** CLEAN UP PANEL DATA SET ***;

* Store labels;
preserve;
collapse (first) var_lab, by(var_str);
local numlab = `=_N';
forvalues j=1/`numlab' {;
	local var_`j' = "`=var_str[`j']'";
	local lab_`j' = "`=var_lab[`j']'";
};
restore;
drop var_lab;

* Reshape to panel;
rename value v_;
reshape wide v_, i(geo geo_str time) j(var_str) string;
rename v_* *;

* Apply labels to variables;
forvalues j=1/`numlab' {;
	label var `var_`j'' "`lab_`j''";
	if length("`lab_`j''")>80 {;
		note `var_`j'': `lab_`j''; // Store note if label is too long;
	};
};

* Save;
xtset geo time;
order geo time;
compress;
save output/global, replace;

