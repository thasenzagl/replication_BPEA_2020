# delimit ;
*clear all;
*macro drop _all;
cap log close;
cap file close _all;
cap program drop _all;
set scrollbufsize 1000000;
set more off;


* Summarize cross-country estimation results;
* MPM 2020-01-21;


*** SETTINGS ***;

local model = "`1'"; // Either "condhet" or "skewt";

* Countries to drop;
local drop_geo_condhet = ""; 		// Cond. het. model;
local drop_geo_skewt = "ESP JPN"; 	// Skew-t model;

*** END SETTINGS ***;


*** LOAD MCMC RESULTS ***;

import delim output/fig_tab/global/var/`model'/results.csv, varnames(1);

drop if strpos("`drop_geo_`model''", geo)>0; // Drop desired countries;
drop if var == "const"; // Drop intercept;


*** PLOTS ***;

* Parameters to report;
local params = "mu sigma";
if "`model'" == "skew_t" {;
	local params = "`params' alpha";
};

preserve;

drop if var == "ylag"; // Drop ylag coef for plots;

* Label "outlier" results with country code;
gen outlier_med_mu = geo if abs(med_mu)>0.1;
gen outlier_probpos_mu = "";
gen outlier_probverypos_mu = "";
gen outlier_probveryneg_mu = "";
gen outlier_med_sigma = geo if abs(med_sigma)>0.1;
gen outlier_probpos_sigma = geo if abs(probpos_sigma-0.5)>0.3;
gen outlier_probverypos_sigma = geo if probverypos_sigma>0.4;
gen outlier_probveryneg_sigma = geo if probveryneg_sigma>0.4;
capture {;
	gen outlier_med_alpha = geo if abs(med_alpha)>0.01;
	gen outlier_probpos_alpha = geo if abs(probpos_alpha-0.5)>0.4;
	gen outlier_probverypos_alpha = geo if probverypos_alpha>0.1;
	gen outlier_probveryneg_alpha = geo if probveryneg_alpha>0.1;
};

encode var, gen(var_id); // Numeric values for variables;
su var_id;
local id_max = `r(max)'; // Number of variables;

foreach s in med probpos probverypos probveryneg {; // For each posterior summary measure...;
	
	foreach p of local params {; // For each type of parameter...;
		
		scatter var_id `s'_`p', ylabel(1/`id_max', ang(h) valuelabel labsize(vsmall)) ysc(reverse)
												msize(vsmall) mlabel(outlier_`s'_`p') mlabcolor(maroon) mlabsize(vsmall)
												graphregion(color(white)) ytitle("") xtitle("");
		graph export output/fig_tab/global/var/`model'/plot_`s'_`p'.png, replace;
		graph export output/fig_tab/global/var/`model'/plot_`s'_`p'.eps, replace;
	
	};
	
};

restore;


*** TABLES ***;

* Summary stats of individual coefficients;
foreach p of local params {; // For each type of parameter...;
	
	gen signif_`p' = qlo_`p'>0 | qhi_`p'<0; // Significance indicator;

	* Latex table with summary stats;
	tabout var using output/fig_tab/global/var/`model'/table_`p'.tex,
		sum
		c(N signif_`p' mean med_`p' mean signif_`p' mean probverypos_`p' mean probveryneg_`p')
		f(0 4 2 2 2)
		clab("N" "med" "signif" "P>delta" "P<-delta")
		style(tex)
		replace;
				
};

if "`model'" == "skew_t" {;

	preserve;
	
	* Summary stats of alpha and TVD;
	collapse (first) mean_alphaidx mean_tvd* *_nu, by(geo);
	tabout geo using output/fig_tab/global/var/`model'/table_tvd.tex,
		sum
		c(mean mean_alphaidx mean mean_tvdavg mean mean_tvdstd mean qlo_nu mean med_nu mean qhi_nu)
		f(3 3 3 1 1 1)
		clab("avg(alpha)" "avg(TVD)" "std(TVD)" "Q1(nu)" "med(nu)" "Q3(nu)")
		style(tex)
		replace;
	
	restore;

};

