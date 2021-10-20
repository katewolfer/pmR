# pmR

A package of tools I currently use for LC-MS(/MS data preprocessing, extraction and visualisation. It uses the most current version of XCMS[https://github.com/sneumann/xcms], which after some experimentation is finally working more or less as it was before the package overhaul. If you're still running on any of the versions before XCMS 3 and are struggling to implement it, please get in touch!
  
## Some of the functions in pmR:
  
### pmRsampleRandomizer
Sample list randomization using batches or total randomization (pseudorandom)
  
### pmRfileConvert
Automated conversion of Thermo/Waters rawfiles to .mzML format
With options for tandem MS data extraction
  
### pmRplotTICs
Plotting of TICs (slow - uses xcms functions currently)
Future feature - BPI chromatogram plotting
  
### pmRparameterization
Automated parameter selection for xcms processing of datasets
  
### pmRexpandList
Expands on simple chemical formula lists to obtain monoisotopic masses of common adducts, to match masses in processed LC-MS data  
  
### pmRmsN
Exploration of MSn data with user input for optimal data extraction  
Database matching - AADB, Pubchem, LipidMaps, HMDB, KEGG  
  
### general statistics
unvariate statistics in conjunction with metadata:  
boxplots, ROC curves, Cohen's D, cross-validation, correlation, linear regression  
PCA  
PLS  
multiple comparisions and FDR  
multiple linear regression  
multiblock analysis  
correlation heatmaps  
Volcano plots  
metadata statistics  
checking normality of feature distributions for univariate testing  
  
## To add:
- van Krevelen plots  
- double bond equivalent plots, H/C and O/C ratio plots
