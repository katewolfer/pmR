pmRMSnExplorer <- function(wd, QCindices) {
  
  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################
  
  ## from xcms vignette 
  ## https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcmsMSn.html
  
  # select files to examine
  fileList <- list.files(wd)
  fileSelection <- seq(from = 1, to = length(list.files(wd)), by = 5)
  fileSelection <- c(fileSelection, length(list.files(wd)))
  
  ## import file using xcmsRaw
  checkFile <- xcmsRaw(fileList[9], profstep = 0.1, includeMSn=TRUE)
  peaks <- findPeaks(checkFile, method="MS1")
  plotTree(xfrag, xcmsFragmentPeakID = 6, textOnly = TRUE)
  
  xs <- xcmsSet(fileList[c(8:11)], method = "MS1")
  xfrag <- xcmsFragments(xs)
  plotTree(xfrag, xcmsFragmentPeakID = 10, textOnly = TRUE)
  
}