pmRparameterization <- function(wd, QCindices) {
  
  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################
  
  # select files to examine - every 5th file here - ensure final file is included
  fileList <- list.files(wd)
  fileSelection <- seq(from = 1, to = length(list.files(wd)), by = 5)
  fileSelection <- c(fileSelection, length(list.files(wd)))
  
  ## import file using xcmsRaw
  checkFile <- xcmsRaw(fileList[9], profstep = 0.1, includeMSn=FALSE)
  
  # get number of cores for parallelization
  nCores = round(detectCores(all.tests = FALSE, logical = TRUE)*0.75)
  BPPARAM=SnowParam(nCores)  

  ## get minimum intensity value for S/N calculation
  # getMinIntensity <- min(checkFile@tic)
  # if (getMinIntensity == 0){
  #   getMinIntensity <- summary(checkFile@tic)[2]
  # }
  
  intensity <- floor(as.numeric(unlist(summary(checkFile@tic)[2])))
  
  # scan numbers for prefilter
  scans <- ceiling(max(checkFile@acquisitionNum)/max(checkFile@scantime))
  if (scans < 5){
    scans <- 5
  }
  
  # minimum chromatographic peak width in s
  CPWmin <- round(scans/(max(checkFile@acquisitionNum)/max(checkFile@scantime)),
                  2)
  
  # maximum chromatographic peak width in s - NEEDS CALCULATION
  CPWmax <- 25
  
  ##
  checkParams <- xcmsSet(#fileList[fileSelection],
                         method = "centWave", 
                         peakwidth = c(CPWmin,CPWmax), 
                         ppm = 10, 
                         integrate = 2, 
                         BPPARAM = BPPARAM,
                         mzCenterFun="wMean",
                         mzdiff = -0.001)
  
  showError(checkParams)
  
  ## check this later - profstep for obiwarp, other settings
  ## Default gap penalties: (gapInit, gapExtend) [by distFunc type]: 
  ## 'cor' = '0.3,2.4' 'cov' = '0,11.7' 'prd' = '0,7.8' 'euc' = '0.9,1.8' 
  start_time <- Sys.time()
  obagds <- retcor(checkParams, 
                   method="obiwarp", 
                   profStep=0.1, 
                   plottype="deviation", 
                   center=5)
  end_time <- Sys.time()
  end_time - start_time
  
  start_time <- Sys.time()
  obagds_01 <- retcor(checkParams, 
                   method="obiwarp", 
                   profStep=0.01, 
                   plottype="deviation", 
                   center=5)
  end_time <- Sys.time()
  end_time - start_time
  
  #  Group peaks together across samples
  bw = 6 # adjust
  MZerror = 0.001

  gagds6 <- group(obagds, method="density", mzwid=MZerror, bw=6) 
  gagds10 <- group(obagds, method="density", mzwid=MZerror, bw=10) 
  gagds3 <- group(obagds, method="density", mzwid=MZerror, bw=3) 
  
  length(intersect(gagds3@groups[,4], gagds10@groups[,4]))
  length(intersect(gagds3@groups[,1], gagds10@groups[,1]))
  
  maxRTshift <- max(gagds@groups[,6] - gagds@groups[,5])
                         
  # # retention time shift in s
  # RTerror = 0
  # 
  # # m/z movement in Da
  # MZerror = 0
    
  # finalised parameters and output
  output <- list('CPWmin' = CPWmin,
                 'CPWmax' = CPWmax,
                 'ppm' = ppm,
                 'LM' = F,
                 'integrate' = 2,
                 'BPPARAM'= BPPARAM,
                 'mzCenterFun' = "wMean",
                 'prefilter' = c(scans, intensity))
  
  return(output)
  
}

# S/N - take beginning/end/other segments and average the baseline
# then 

# try response surface model for parameters selected on random subset?
# then visualize results?
# model the baseline?

# removal of background/contaminant peaks (stuff that has low CV% across ALL
# samples? )
