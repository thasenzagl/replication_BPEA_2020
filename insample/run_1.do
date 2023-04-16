clear all
set more off


* Replication files: When is Growth at Risk?
* Step 1: Create data and benchmark results


* Create data sets used for further analysis
cd data/global
do create_data
cd ../us
do create_data

* Create crisis statistics cited in Section IV.A
cd ../global
do crises

* Table S.4
cd ../../estim
do benchmark_regr

cd ..
