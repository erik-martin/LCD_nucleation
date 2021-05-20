Each NaCl concentration directory contains the file:
data.mat

This file contains the averaged saxs data as a function of time in Matlab matrix format. When read, the file is organzed as a cell array where each cell contains all the measurements for each time point. The cells are organized as 1211x2xX matrices where the 3 columns are q,I(q), and error. X represents tthe number of individual measurements in the experiment. 

The .mat files can be read by the analysis scripts:
	Analyze_laminarflow_data.m
	plot_raw_data.m

fuzzy_colloid.m - generates the simulated data using the fuzzy colloid model
mixing.csv - contains the mixing time simulations
