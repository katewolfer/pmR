pmRMSnExplorer <- function(wd, QCindices) {
  
  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################
  
  ## from xcms vignette 
  ## https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcmsMSn.html
  
  # select files to examine
  wd <- getwd()
  fileList <- list.files(wd)
  fileSelection <- seq(from = 1, to = length(list.files(wd)), by = 5)
  fileSelection <- c(fileSelection, length(list.files(wd)))
  
  ## import file using xcmsRaw
  checkFile <- xcmsRaw(fileList[10], profstep = 0.1, includeMSn=TRUE)
  peaks <- findPeaks(checkFile, method="MS1")
  
  xs <- xcmsSet(fileList[c(10:13)], method = "MS1")
  xfrag <- xcmsFragments(xs)
  
  plotTree(xfrag, xcmsFragmentPeakID = 6, textOnly = TRUE)
  plotTree(xfrag, xcmsFragmentPeakID = 10, textOnly = TRUE)
  
}