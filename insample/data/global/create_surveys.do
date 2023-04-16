# delimit ;


* Create consumer/business survey data set;
* MPM 2020-01-16;


*** SETTINGS ***;

* Files;
local xls_file = "surveys/Survey data 10_01_2020.xlsx"; // Spreadsheet with survey data;
local xls_sheets = "Australia Belgium Canada Switzerland Germany Spain France UK Italy Japan Netherlands Sweden US"; // Sheets to load from in spreadsheet;

* See also additional settings under "VARIABLE CHOICES" below;

*** END SETTINGS ***;


*** LOAD MONTHLY SERIES ***;

cap erase "output/surveys.dta";

foreach c of local xls_sheets {;

	import excel using "`xls_file'", sheet("`c'") clear;
	
	* Drop empty variables;
	foreach v of varlist * {;
		if "`=`v'[3]'" == "." {;
			drop `v';
		};
	};
	
	* Go through remaining variables in pairs (time and value);
	ds;
	local vars = "`r(varlist)'";
	local numvars: word count `vars';
	
	forvalues j=2(2)`numvars' {; // For every second variable in memory...;
		
		local the_series_code: word `j' of `vars'; // Stata variable name (values);
		local the_prev_series_code: word `=`j'-1' of `vars'; // Stata variable name of previous variable (time variable);
		local the_series_name = "`=`the_series_code'[2]' `=`the_series_code'[1]'"; // Real name of series;
		
		preserve;
		
		* Keep only time and values for the current series;
		rename `the_prev_series_code' time_str;
		rename `the_series_code' value;
		keep time_str value;
		drop if time_str=="";
		
		* Series and country names;
		gen series = "`the_series_name'";
		gen country = "`c'";
		
		* Append to existing data;
		cap append using output/surveys;
		cap save output/surveys, replace;
		while _rc != 0 {;
			sleep 100;
			cap save output/surveys, replace;
		};
		
		restore;
		
	};

};

use output/surveys, clear;

* Time variable;
gen time_m = monthly(time_str, "YM");
format time_m %tm;
drop time_str;

* Values to numeric;
destring value, force replace;
drop if missing(value);

* Country names to ISO code;
kountry country, from(other) stuck;
rename _ISO3N_ country2;
kountry country2, from(iso3n) to(iso3c);
rename _ISO3C_ geo_str;
drop country*;

* Import country IDs;
save output/surveys, replace;
use output/oecd_eo, clear;
collapse (first) geo, by(geo_str);
merge 1:m geo_str using output/surveys, nogen keep(match); // Load data again;
label val geo geo;


*** VARIABLE CHOICES ***;

gen var_str = "";

* For each country, manually select variables of the types: economic/business sentiment, consumer sentiment, PMI;

* AUS;
replace var_str = "CONSSENT" if geo_str == "AUS" & series == "Consumer Sentiment Westpac-Melbourne Institute";
replace var_str = "PMI" if geo_str == "AUS" & series == "Mfg PMI AIG/PWC";

* BEL;
replace var_str = "CONSSENT" if geo_str == "BEL" & series == "Consumer Confidence EC";
replace var_str = "ECONSENT" if geo_str == "BEL" & series == "Busines Survey: Mfg BNB";

* CHE;
replace var_str = "PMI" if geo_str == "CHE" & series == "PMI Procure.ch";

* ESP;
replace var_str = "CONSSENT" if geo_str == "ESP" & series == "Consumer Confidence EC";

* FRA;
replace var_str = "ECONSENT" if geo_str == "FRA" & series == "Business Climate INSEE";

* GBR;
replace var_str = "ECONSENT" if geo_str == "GBR" & series == "Economic Sentiment Indicator EC";
replace var_str = "CONSSENT" if geo_str == "GBR" & series == "Consumer Confidence GfK";

* ITA;
replace var_str = "CONSSENT" if geo_str == "ITA" & series == "Consumer Confidence ISTAT";

* JPN;
replace var_str = "CONSSENT" if geo_str == "JPN" & series == "Consumer Confidence Cabinet Office";

* NLD;
replace var_str = "ECONSENT" if geo_str == "NLD" & series == "Economic Sentiment EC";
replace var_str = "CONSSENT" if geo_str == "NLD" & series == "Consumer Confidence EC";

* SWE;
replace var_str = "ECONSENT" if geo_str == "SWE" & series == "Economic Sentiment EC";
replace var_str = "CONSSENT" if geo_str == "SWE" & series == "Consumer Confidence NIER";
replace var_str = "PMI" if geo_str == "SWE" & series == "Mfg PMI Swedbank";

* USA;
replace var_str = "ECONSENT" if geo_str == "USA" & series == "Business Outllook Survey Philly Fed";
replace var_str = "CONSSENT" if geo_str == "USA" & series == "Consumer Confidence Conference Board";
replace var_str = "PMI" if geo_str == "USA" & series == "Mfg PMI ISM";

drop if var_str=="";


*** COLLAPSE TO QUARTERLY AVERAGES ***;

gen time = qofd(dofm(time_m));
format time %tq;
collapse (mean) value, by(geo geo_str time var_str);


*** RESHAPE ***;

* Reshape to panel;
reshape wide value, i(time geo_str geo) j(var_str) string;
rename value* *;

* Set panel;
xtset geo time;


*** LABEL AND SAVE ***;

label var ECONSENT "Economic/business sentiment/confidence";
label var CONSSENT "Consumer sentiment/confidence";
label var PMI "Purchasing Managers Index";

order geo geo_str time;
compress;
save output/surveys, replace;


*** MERGE WITH OTHER DATA ***;

merge 1:1 geo time using output/global, nogen;
xtset geo time;

* Save;
compress;
save output/global, replace;



