pmRparameterization <- function(keyFile, QCindices) {

  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################

  # select files to examine
  fileList <- list.files(getwd(), pattern = "cdf")

  ## import file using xcmsRaw
  checkFile <- xcmsRaw(fileList[3], profstep = 0.1, includeMSn=FALSE)

  ######################
  ## Plotting of TICs ##
  ######################

  # ## get file list and types
  # startList <- 1
  # endList <- 56
  # keyFile <-  read.csv("eliquid experiments May 2020.csv")
  # trimFilenames <- paste(gsub( ".mzML", "", fileList[c(startList:endList)]))
  # findFiles <- which(is.element(keyFile$fileName, trimFilenames))
  # keyFile$sample[findFiles]

  # slow - uses xcms functions currently
  # Future feature - BPI chromatogram plotting
  source("pmRplotTICs.R")
  TICplot <- pmRplotTICs(fileList[2])
  backupPlot <- TICplot # back up the original plot
  TICplot

  # reduce windows for viewing
  TICplot <- TICplot + coord_cartesian(ylim = c(0,1.5e9))
  TICplot <- TICplot + coord_cartesian(xlim = c(0,3))

  ggsave("TICplot_14Sept2020.png", plot = TICplot, device = NULL, path = NULL,
         scale = 1, width = 32, height = 16, dpi = 600,units = c("cm"))

  # for single file
  ## produce TIC plot
  pullData <- as.data.frame(cbind(checkFile@scantime, checkFile@tic))
  colnames(pullData) <- c("RT", "intensity")

  p1 <- ggplot(pullData, aes(x=RT, y=intensity))
  p1 <- p1 + geom_line(na.rm = TRUE)
  p1 <- p1 +  xlab('RT (s)') + ylab('intensity')
  p1 <- p1 + theme(plot.title = element_text(hjust = 0.5, size = 18))
  p1 <- p1 + theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour="light grey", size=0.01),
                   panel.border = element_blank(),
                   panel.grid.major = element_line(colour="light grey", size=0.01),
                   panel.background = element_blank())
  p1

  # pick peaks from data for assessment
  peak <- findPeaks.matchedFilter(checkFile,step = 0.02, mzdiff = 0.001) # scanrange=c(50,70))
  #peak <- findPeaks(checkFile,method = "centWave")
  #peak <- findPeaks.matchedFilter(checkFile,scanrange=c(100,400)) # for profile data
  #plotPeaks(checkFile, peak)
  peakData <- as.data.frame(peak@.Data)
  peakData <- peakData[order(peakData$rt),]

  BPI <- peakData[,c("rt", "into")]
  BPI$duplic <- duplicated(BPI$rt)
  BPI$remove <- FALSE

  for (i in 1:nrow(BPI)){
    if (BPI$duplic[i] == FALSE){
    } else {
      findAll <- which(BPI$rt == BPI$rt[i])
      findMax <- max(BPI$into[findAll])
      findKeeper <- which(BPI$into[findAll] == findMax)
      BPI$remove[findAll[-findKeeper]] <- TRUE
    }
  }

  BPI <- BPI[-which(BPI$remove == TRUE),]
  BPI$rt <- BPI$rt/60
  plot(BPI$rt, BPI$into)

  # ## plot the difference between into and intb in the peak picking
  # p <- ggplot(data = peakData) + geom_line(data = peakData,
  #                                          aes(x = as.numeric(rownames(peakData)),
  #                                          y = peakData[,7]), colour = "blue")
  # p <- p + geom_line(data = peakData,
  #                   aes(x = as.numeric(rownames(peakData)),
  #                       y = peakData[,8]), colour = "orange")
  #
  # p <- p + theme(axis.line = element_line(colour = "black"),
  #                  panel.grid.minor = element_line(colour="light grey", size=0.01),
  #                  panel.border = element_blank(),
  #                  panel.grid.major = element_line(colour="light grey", size=0.01),
  #                  panel.background = element_blank())
  # p <- p + xlab("Feature index (arbitrary)")
  # p <- p + ylab("Intensity")
  # p

  # ## try to make a BPI
  # peakInfo <- peakData[c("rt", "into")]
  # peakInfo <- peakInfo[order(peakInfo$rt),]
  # plot(peakInfo$rt, peakInfo$into, type = 'l')
  #
  # # https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
  # raw_data <- readMSData(files = fileList, pdata = new("NAnnotatedDataFrame", pd),
  #                        mode = "onDisk")
  #
  # bpis <- chromatogram(checkFile, aggregationFun = "max")
  #
  # # get number of cores for parallelization
  # totalCores <- detectCores(all.tests = FALSE, logical = TRUE)
  # nCores <- round(totalCores*0.75)
  # BPPARAM <- SnowParam(nCores)

## Signal to noise calculation

  # getMinIntensity <- min(checkFile@tic)
  # if (getMinIntensity == 0){
  #   getMinIntensity <- summary(checkFile@tic)[2]
  # }

  # min and max signals, relate to mean/median signal intensity
  baseSignal <- floor(as.numeric(unlist(summary(checkFile@tic)[2])))
  maxSignal <- ceiling(as.numeric(unlist(summary(checkFile@tic)[6])))

  plot(as.numeric(paste(unlist(summary(checkFile@tic)))))

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

  ## actually start processing whole dataset here - select peak picking
  ## depending on centroid or profile data

  ## CENTROIDED
  # checkParams <- xcmsSet(fileList,
  #                        #method = "centWave",
  #                        method = "matchedFilter",
  #                        peakwidth = c(CPWmin,CPWmax),
  #                        ppm = 2,
  #                        integrate = 2,
  #                        BPPARAM = BPPARAM,
  #                        mzCenterFun="wMean",
  #                        mzdiff = -0.0001)

  ## PROFILE
  checkParams <- xcmsSet(fileList,
                         method = "matchedFilter",
                         BPPARAM = BPPARAM,
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
                   center=1)
  end_time <- Sys.time()
  end_time - start_time

  # start_time <- Sys.time()
  # obagds_01 <- retcor(checkParams,
  #                  method="obiwarp",
  #                  profStep=0.01,
  #                  plottype="deviation",
  #                  center=5)
  # end_time <- Sys.time()
  # end_time - start_time

  #  Group peaks together across samples
  bw = 10 # adjust
  MZerror = 0.001
  gagds <- group(obagds, method="density", mzwid=MZerror, bw=bw)

  # length(intersect(gagds3@groups[,4], gagds10@groups[,4]))
  # length(intersect(gagds3@groups[,1], gagds10@groups[,1]))
  # maxRTshift <- max(gagds@groups[,6] - gagds@groups[,5])

  # replace peaks with zero intensity with background value
  fds <- fillPeaks(gagds)
  fill <- fds
  vals.raw <- groupval(fill, "medret", "into")
  vals.raw.meta <- cbind(groups(fill), vals.raw)
  write.csv(vals.raw.meta, file="extracted_peaks.csv")

  # # finalised parameters and output
  # output <- list('CPWmin' = CPWmin,
  #                'CPWmax' = CPWmax,
  #                'ppm' = ppm,
  #                'LM' = F,
  #                'integrate' = 2,
  #                'BPPARAM'= BPPARAM,
  #                'mzCenterFun' = "wMean",
  #                'prefilter' = c(scans, intensity))
  #
  # return(output)
  #


  ## View annotated peaks

  annotatedPeaks <- read.csv("eliquid annotated peaks.csv")
  annotatedPeaks <- annotatedPeaks[which(annotatedPeaks$putative.ID != ""),]

  rtr <- c(annotatedPeaks$rtmin[i], annotatedPeaks$rtmin[i])
  mzr <- c(annotatedPeaks$mzmin[i], annotatedPeaks$mzmin[i])
  plot.xcmsEIC(fds, rtrange = rtr)
  xeic.corrected <- getEIC(fds, rt = "corrected", groupidx = 1:nrow(fds@groups))

  i = 5
  plot(xeic.corrected, fds,
       groupidx=which(paste(rownames(vals.raw.meta)) == annotatedPeaks$tag[i]),
       col = c("orange", "blue"), main = NULL)

  title("change this")

}

# S/N - take beginning/end/other segments and average the baseline
# then

# try response surface model for parameters selected on random subset?
# then visualize results?
# model the baseline?

# removal of background/contaminant peaks (stuff that has low CV% across ALL
# samples? )
