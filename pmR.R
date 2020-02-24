#####################################
## pmR package                     ##
## Kate Wolfer, Universitaet Basel ##
## v1.0, 03 October 2019           ##
#####################################

## check packages are installed and install if not
## need to set the library location so there aren't issues later on
require(ggplot2)
require(xcms)
require(BiocParallel)
require(RColorBrewer)
require(reshape2)

## for installing BiocParallel and XCMS
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("xcms", lib="C:/Program Files/R/R-3.5.3/library")
# BiocManager::install("BiocParallel", lib="C:/Program Files/R/R-3.5.3/library")


## set Proteowizard program path
#proteowizardPath <- 'C:/Program Files/ProteoWizard/ProteoWizard 3.0.19276.7124c5404'
proteowizardPath <- "C:/Users/wolfer0000/AppData/Local/Apps/ProteoWizard 3.0.20051.87e974e37"# 64-bit"

proteowizardPath <- normalizePath(file.path('C:',
                                            'Users',
                                            'wolfer0000',
                                            'AppData',
                                            'Local',
                                            'Apps',
                                            #'ProteoWizard 3.0.20051.87e974e37'))
                                            'ProteoWizard 3.0.20051.87e974e37 64-bit'))


## set target folder where raw files live
targetFolder     <- normalizePath(file.path('C:',
                                            'Users',
                                            'wolfer0000',
                                            'Desktop',
                                            'test_data',
                                            'baseline'))

# create new directory for mzML converted data
dir.create(paste0(targetFolder, '\\converted_data'))
saveFolder        <- paste0(targetFolder, '\\converted_data')

# perform file conversion
source("pmRfileConvert.R")
pmRfileConvert(proteowizardPath, targetFolder, saveFolder)

# set the working directory to the location of the converted files
setwd(saveFolder)

## plot all TICs
fileList <- list.files(getwd())
TICplot <- pmRplotTICs(fileList)

# parameters
params <- pmRparameterization(getwd(), 10)


# top50 <- sort(checkFile@tic, decreasing = TRUE)
# summary(top50[1:(round(length(checkFile@tic)/2))])
# summary(top50[round(length(checkFile@tic)/2):length(top50)])

