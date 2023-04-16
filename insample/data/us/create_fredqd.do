# delimit ;


* Create FRED-QD data set;
* MPM 2020-01-20;


*** SETTINGS ***;

* Files;
local file_fredqd = "fredqd/fredqd_20191216.csv"; // FRED-QD data;
local file_fredqdser = "fredqd/fredqd_series.xlsx"; // Selected FRED-QD series;

*** END SETTINGS ***;


*** LOAD DATA ***;

* Import data from CSV;
import delim using `file_fredqd', varnames(1) case(preserve);
rename * val_*;
rename val_sasdate date;
reshape long val_, i(date) j(fred_mnemonic) string;
save output/fredqd, replace;

* Import series;
import excel using `file_fredqdser', sheet("series") firstrow clear;
keep if use==1;
merge 1:m fred_mnemonic using output/fredqd, nogen keep(match);


*** CLEAN UP ***;

rename val_ val;
drop if missing(val);
drop use;
replace sw_id = "" if sw_id == "n.a.";

gen time = qofd(date(date, "MDY"));
format time %tq;
drop if missing(time);
drop date;

rename (fred_mnemonic category) (fred_mnemonic_str category_str);

destring, replace;
xtset id time;


*** MANUAL VARIABLE DEFINITIONS ***;

* Interest rates used to calculate spreads;
foreach v in FEDFUNDS TB3MS GS10 {;
	gen `v'_val = val if fred_mnemonic_str == "`v'";
	egen `v' = max(`v'_val), by(time);
};

* Transform TB3MS to spread;
replace val = val - FEDFUNDS if fred_mnemonic_str == "TB3MS";
replace tc = 1 if fred_mnemonic_str == "TB3MS";
replace fred_mnemonic_str = "TB3MSFEDFUNDSx" if fred_mnemonic_str == "TB3MS";

* Transform GS10 to spread;
replace val = val - TB3MS if fred_mnemonic_str == "GS10";
replace tc = 1 if fred_mnemonic_str == "GS10";
replace fred_mnemonic_str = "GS10TB3MSx" if fred_mnemonic_str == "GS10";

* Transform AAA & BAA to spreads;
foreach v in AAA BAA {;
	replace val = val - GS10 if fred_mnemonic_str == "`v'";
	replace tc = 1 if fred_mnemonic_str == "`v'";
	replace fred_mnemonic_str = "`v'GS10x" if fred_mnemonic_str == "`v'";
};

drop FEDFUNDS TB3MS GS10;


*** TRANSFORM TO STATIONARITY ***;

gen valt = val if tc == 1;
replace valt = D.val if tc == 2;
replace valt = D2.val if tc == 3;
replace valt = 100*log(val) if tc == 4;
replace valt = 100*(log(val)-log(L.val)) if tc == 5;
replace valt = 100*((log(val)-log(L.val))-(log(L.val)-log(L2.val))) if tc == 6;
replace valt = 100*(val/L.val-1) if tc == 7;

* Save;
encode fred_mnemonic_str, gen(fred_mnemonic);
encode category_str, gen(category);
order fred_mnemonic category time;
compress;
save output/fredqd, replace;

