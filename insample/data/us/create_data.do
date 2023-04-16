# delimit ;
clear all;
macro drop _all;
cap log close;
cap file close _all;
cap program drop _all;
set scrollbufsize 1000000;
set more off;
log using output/create_data.log, replace;


* Create U.S. data set;
* MPM 2019-12-17;


*** SETTINGS ***;

* Files;
local dat_global = "../global/output/global.dta"; // Global database (for loading certain extra U.S. variables);
local file_export = "output/data_export.csv"; // Export file;

* Variables to import from global database;
local global_vars = "COMMCRB ECONSENT CONSSENT PMI STOCKVOL";
local global_vars_dlog = "COMMCRB STOCKVOL"; // Transform these to log differences;

* Export sample;
local export_sample = "1975q2,2019q2";

*** END SETTINGS ***;


*** CREATE FRED-QD DATA ***;

do create_fredqd;


*** COMBINE WITH U.S. DATA FROM GLOBAL DATABASE ***;

* Reshape FRED-QD data;
use output/fredqd;
keep time fred_mnemonic_str valt;
reshape wide valt, i(time) j(fred_mnemonic_str) string;
rename valt* *;

* Load data from global database;
gen geo_str = "USA";
merge 1:1 time geo_str using `dat_global', nogen keepusing(`global_vars') keep(match master);
drop geo_str;
tsset time;

* Transform to log growth rates;
foreach v of varlist `global_vars_dlog' {;
	gen `v'_dlog = log(`v'/L.`v');
	drop `v';
	rename `v'_dlog `v';
};


*** RENAME VARIABLES ***;

rename ?AAGS10x ?AASPR;
rename AHETPIx EARNINGS;
rename AMDMNOx ORDERNEW;
rename AMDMUOx ORDERUNFIL;
rename BUSLOANSx LOANSCORP;
rename CONSUMERx LOANSHH;
rename CPF3MTB3Mx CPAPERSPR;
rename DPIC96 DISPINC;
rename EXPGSC1 EXPORT;
rename GCEC1 CONSGOVT;
rename GDPC1 GDP;
rename GFDEGDQ188S DEBTGOVT;
rename GPDIC1 INVESTM;
rename GS10TB3MSx TERMSPR;
rename HOANBS HOURS;
rename HOUST HOUSESTART;
rename IMPGSC1 IMPORT;
rename INVCQRMTSPL INVENTO;
rename PAYEMS EMPL;
rename PCECC96 CONSPRIV;
rename PCEPILFE PCEPRICE;
rename PERMIT HOUSEPERMIT;
rename RSAFSx RETAIL;
rename SP500 STOCKPRICE;
rename SPdivyield DIVYIELD;
rename TB3MSFEDFUNDSx SHORTSPR;
rename TCU CAPUTIL;
rename TLBSHNOx LIABHH;
rename TLBSNNCBx LIABCORP;
rename TNWBSHNOx NWHH;
rename TNWMVBSNNCBx NWCORP;
rename TWEXMMTH EXCHTRW;
rename ULCNFB ULC;
rename USSTHPI HOUSEPRICE;
rename VXOCLSx VXO;

* Save combined data;
keep if tin(`export_sample');
compress;
save output/data_export, replace;


*** EXPORT ***;

gen date = dofq(time);
format date %td;
order *, alpha;
order date;
xtset, clear;
drop time;
export delim `file_export', replace;

use output/data_export, clear;


log close;
