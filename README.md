# Detzel--Novy-Marx--and-Velikov
 Code used to create results in Detzel, Novy-Marx, and Velikov (JF, Forthcoming), Model Comparison with Transaction Costs

This repository contains code used to create the results in Detzel, Novy-Marx, and Velikov (JF, Forthcoming), Model Comparison with Transaction Costs. This code is to be used in conjunction with the MATLAB asset pricing package that accompanies Novy-Marx and Velikov (WP, 2022), Assaying Anomalies. 

The order of operations to replicate the results in Detzel, Novy-Marx, and Velikov (WP, 2022) is:

1. Download and follow the instructions for setting up the MATLAB asset pricing package from https://github.com/velikov-mihail/AnomalyCookbookOfficial.git
	* The results in Detzel, Novy-Marx, and Velikov (JF, Forthcoming) use the pre-release v0.2 of the MATLAB asset pricing package, which will made public soon. The code is also available from the authors upon request.
3. Download the code in this repository.
4. Run dnmv.m. The script requires setting up the directories for the MATLAB asset pricing package repository and this repository. It starts a log file and calls multiple other scripts which perform the following functions:  
	* make_main_factors_tcosts.m downloads the factors from Hou, Xue, and Zhang (RFS, 2015) and Asness and Frazzini (JPM, 2013), replicates all factors, calculates the transaction costs for, and stores all factors  
	* make_mitigated_factors.m creates and stores mitigated versions of all factors  
	* make_freak_factor.m creates and stores the "freak" factor based on low-volatility industry relative reversals used in Figure 1  
	* organize_factors.m loads all factors, organizes them in a tidy MATLAB structure, and stores it  
	* run_full_sample_results.m creates and stores a structure with results for the full-sample results
	* run_bootstraps.m creates and stores structures with results for the bootstrap results
 	* run_oos_results.m creates and stores structures with results for the out-of-sample MVE portfolio results
	* run_anomalies_results.m creates and stores structures with results for the anomaly section
	* make_tables.m prints all tables in ~/Results/DetzelNovy-MarxVelikovTablesOutput.txt
	* make_figures.m creates and stores all figures in ~/Figures/
   
