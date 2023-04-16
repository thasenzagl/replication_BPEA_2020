clear all
set more off


* Replication files: When is Growth at Risk?
* Step 4: Create extra figures and tables


cd estim

* Figures 13, S.14 and Tables 2-3
do summ_results_global condhet

* Tables 4, S.5
clear all
do summ_results_global skewt

cd ..