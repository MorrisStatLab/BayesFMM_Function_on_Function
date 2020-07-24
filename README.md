# FFR
Function-on-Function Regression
Matlab files for 
'Function-on-Function Regression for the Identification of Epigenetic Regions Exhibiting Windows of Susceptibility to Environmental Exposures' by Zemplenyi, Meyer, Cardenas, et al. (pre-print available on arXiv: https://arxiv.org/abs/1912.07359.
This code base builds upon the WFMM package created at MD Anderson, available here: https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/70.

Files:
* run_FFR.m: script to initiate a function-on-function regression analysis. You can use either the simulated exposure ("simX_N400_T90.csv") and response ("simY_N400_S100.csv") data provided or specify the exposure, response, and covariate data for your custom analysis. 
* generate_X.m: script used to generate the simulated exposure data contained in "simX_N400_T90.csv." 
* generate_Y.m: script used to generate the simulated response data contained in "simY_N400_S100.csv." 
* make_sim_heatmaps.m: function called by script "run_FFR.m" that creates three heatmaps (estimated beta surface, Bayesian false discovery rate, and simutaneous band scores) using the "results.mat" object from the FFR analysis. 
* The following outputs are sent to the "Results" folder after the script has finished:
	* "FFR_output.mat" raw results object 
	* three heatmaps (estimated beta surface, Bayesian False Discovery Rate inferential procedure, and Simultaneous Band Score inferential procedure) 

