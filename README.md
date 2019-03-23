
This is repo for code used in generating final report for CSE599, data analysis in system research.

The repo is organized as the following (introduced in chronical order if one wants to reproduce the analysis):
./data_generation_code/
	-  fix_size_logistic.R and random_size_logistic.R generates RData files (which contains the data frame) from the Monte Carlo simulation
	- logistic_analysis_fixed_size.R/logistic_analysis_random_size.R takes the RData files generated in the last step as input and outputs csv files in  the analysis_input folder

./analysis_input/
	- contains csv files for different summary statistics and effect size measures computed using the logistic_analysis_{fixed/random}_size.R codes

./analysis_code/
	- plot_fixed_random_df.R plots every figure correlation related in the report
	- plot_fixed_random_alt_measure.R plots every figure linear regression related in the report
	- plot_fixed_random_alt_measure_or.R plots every figure logistic regression related in the report

./report/
	- contains the figure and tex used in the final report