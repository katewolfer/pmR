pmRsetup <- function(){

  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################

  ## check packages are installed and install if not
  ## need to set the library location so there aren't issues later on
  fetchLibrary = "C:/Program Files/R/R-4.0.2/library"
  
  # minimum required packages
  require(BiocParallel)
  require(BiocManager)
  require(xcms)
  require(MSnbase)

  # packages for plotting
  require(RColorBrewer)
  require(reshape2)
  require(plyr)
  require(ggplot2)
  require(ggrepel)
  require(ggsignif)
  require(gridExtra)

  ## for installing BiocParallel and XCMS
  packages <- c("BiocManager", "xcms", "BiocParallel", "MSnbase", "multtest")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", lib = fetchLibrary)
    BiocManager::install("xcms", lib = fetchLibrary)
    BiocManager::install("BiocParallel", lib = fetchLibrary)
    BiocManager::install("MSnbase", lib = fetchLibrary)
    BiocManager::install("multtest", lib = fetchLibrary)
  }
}
