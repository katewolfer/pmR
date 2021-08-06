## pmR extra code repo

################
## TO DO list ##
################

## pmRmsN
# Exploration of MSn data with user input for optimal data extraction
# database matching - AADB, Pubchem, LipidMaps, HMDB, KEGG

## pmRquant
# Quantification of features
# Checking run order effects

### general statistics
# unvariate statistics in conjunction with metadata:
# boxplots, ROC curves, Cohen's D, cross-validation, correlation, linear regression
# PCA
# PLS
# multiple comparisions and FDR
# multiple linear regression
# multiblock analysis
# correlation heatmaps
# Volcano plots
# metadata statistics
# checking normality of feature distributions for univariate testing

## pmRsampleRandomizer
# Sample list randomization using batches or total randomization (pseudorandom)



####################
## pmRfileConvert ##
####################

convertFiles <- FALSE

if (convertFiles == TRUE){
  ## set Proteowizard program path




  targetFolder     <- normalizePath(file.path('C:',
                                              'Users',
                                              'Kate',
                                              'switchdrive',
                                              'data',
                                              'nicotine_analysis'))

  targetFolder    <- normalizePath(file.path('C:',
                                             'Users',
                                             'Kate',
                                             'switchdrive',
                                             'data',
                                             'nicotine_analysis'))

  # create new directory for mzML converted data
  dir.create(paste0(targetFolder, '\\converted_data'))
  saveFolder <- paste0(targetFolder, '\\converted_data')

  ## pmRfileConvert
  # Automated conversion of Thermo/Waters rawfiles to .mzML format
  # With options for tandem MS data extraction
  source("pmRfileConvert.R")
  pmRfileConvert(proteowizardPath, targetFolder, saveFolder)

  # set the working directory to the location of the converted files
  setwd(saveFolder)
}


#########################
## pmRparameterization ##
#########################
source("pmRXCMSparameters.R")
getParams <- pmRXCMSparameters()
pmRparameterization(wd, QCindices)


###################
## pmRexpandList ##
###################

# Expands on simple chemical formula lists to obtain monoisotopic masses of
# common adducts. Can be used with (function) to match masses in processed
# LC-MS data
source("pmRexpandList.R")
df <- read.csv("eliquid_compounds_14Sept2020.csv", check.names = FALSE,
               stringsAsFactors = FALSE)
polarityMode <- "pos"
ppm <- 5
getList <- pmRexpandList(df, ppm, polarityMode)

vals.raw.meta <- as.data.frame(vals.raw.meta)
#vals.raw.meta$mzmed

checkPPM <- NULL
filterTable <- NULL

for (i in 1:nrow(vals.raw.meta)){
  for (j in 1:nrow(getList)){
    checkPPM <- vals.raw.meta$mzmed[i] > getList$lower_ppm[j] && vals.raw.meta$mzmed[i] < getList$upper_ppm[j]
    if (checkPPM == TRUE) {
      slotRow <- cbind(vals.raw.meta[i,],getList[j,])
      filterTable <- rbind(filterTable, slotRow)
    }
  }
}

write.csv(getList, 'all_eliquid_masses_14Sept2020.csv')



## ------------------------------------------------------------##


