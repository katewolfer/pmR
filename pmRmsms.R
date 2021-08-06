############################################################
## MS/MS extraction, processing and visualisation testing ##
## Kate Wolfer, Universitaet Basel                        ##
## July 2021                                              ##
############################################################

## Useful links for mzR/msnBase package functions/data processing:
# https://lgatto.github.io/2020-02-17-RProt-Berlin/raw-ms-data-mzr-and-msnbase.html
# https://shiny.rstudio.com/tutorial/written-tutorial/lesson2/
# https://lgatto.github.io/MSnbase/articles/v01-MSnbase-demo.html

#rawdata <- readMSData(fileList[1], verbose = FALSE)
#rawdata <- readMSData(fileList[1], msLevel = 2, verbose = FALSE)

## required packages
require(BiocManager)
require(BiocParallel)
require(xcms)
require(ggplot2)
require(shiny)


## get all mzML files in directory
fileList <- list.files(getwd(), pattern = "mzML")

## select the file from the file list to be analysed
selectFile <- 3

## open the file to access the data
rawData <- openMSfile(fileList[selectFile], backend = "pwiz", verbose = FALSE)

## pull out the full experiment summary
checkHeader <- header(rawData)

## this represents all the data channels which have been acquired in the Thermo raw file
getAllPrecursors <- unique(checkHeader$filterString)

## select each of the types of experiment in the file and restructure
getExpt <- 1

## plotting the TIC or BPI of the main channel, if one exists
#if (length(grep("ESI Full ms", getAllPrecursors[getExpt])) > 0){

findExperiment <- checkHeader$acquisitionNum[which(checkHeader$filterString == getAllPrecursors[getExpt])]
exptDetails <- checkHeader[findExperiment,]
exptTitle <- getAllPrecursors[getExpt]

## get TIC
totalIC <- as.data.frame(cbind(exptDetails$seqNum, exptDetails$totIonCurrent))
colnames(totalIC) <- c("rt", "intensity")
totalIC$rt <- as.numeric(totalIC$rt)
totalIC$intensity <- as.numeric(totalIC$intensity)

## get BPI
BPI <- as.data.frame(cbind(exptDetails$seqNum, exptDetails$basePeakIntensity))
colnames(BPI) <- c("rt", "intensity")
BPI$rt <- as.numeric(BPI$rt)
BPI$intensity <- as.numeric(BPI$intensity)

source("pmRplotTIC.R")
chrom <- BPI
getTIC <- pmRplotTIC(chrom)
getFileName <- fileList[selectFile]
getTIC <- getTIC + ggtitle(paste0("Experiment: ", getFileName, "\n", exptTitle))
getTIC

#}

getExpt <- 26

findExperiment <- checkHeader$acquisitionNum[which(checkHeader$filterString == getAllPrecursors[getExpt])]
exptDetails <- checkHeader[findExperiment,]
exptTitle <- getAllPrecursors[getExpt]

collectPeakData <- peaks(rawData, exptDetails[,1])
rl <- rapply(collectPeakData, length, how="list")
rl <- rapply(collectPeakData, max, how="list")

getMasses <- NULL
getRTs <- exptDetails$retentionTime
plot(exptDetails$retentionTime, exptDetails$totIonCurrent, type = 'l')

# lo <- loess(exptDetails$totIonCurrent~exptDetails$retentionTime)
# plot(exptDetails$retentionTime,exptDetails$totIonCurrent)
# lines(predict(lo), col='red', lwd=2)

# https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
smoothingSpline = smooth.spline(exptDetails$retentionTime, log10(exptDetails$totIonCurrent), spar=0.2)
plot(exptDetails$retentionTime, log10(exptDetails$totIonCurrent))
lines(smoothingSpline)

getDensity <- density(exptDetails$basePeakIntensity)
plot(getDensity)

require(Rwave)
cwt(exptDetails$basePeakIntensity, noctave=4, nvoice=50)

for (b in 1:length(collectPeakData)){
  collectPeakData[1][[1]]
}


## spectrum plot
for (getScan in c(1:nrow(exptDetails))){
  chrom <- as.data.frame(peaks(rawData, exptDetails[getScan,1]))
  colnames(chrom) <- c("mz", "intensity")
  chrom$mz <- as.numeric(chrom$mz)
  chrom$intensity <- as.numeric(chrom$intensity)

  source("pmRplotSpectrum.R")
  getSpectrum <- pmRplotSpectrum(chrom)
  getSpectrum <- getSpectrum + ggtitle(paste0("Experiment: ", getFileName, "\n", exptTitle, ", scan ", getScan))
  plot(getSpectrum)
  readline(prompt="Press [enter] to continue")
  dev.off()
}


#ggsave(paste0(fileName,".png"), plot = chromPlot, device = NULL, path = NULL, scale = 1, width = 17, height = 11, dpi = 500,units = c("cm"))


## need to relate scan numbers to actual retention times
# https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm



i <- 247
pi <- peaks(rawData, MS2[i, 1])

## prototyping of peak data extraction and plotting
i <- 247
MS2[i,]
pi <- peaks(rawData, MS2[i, 1])
plot(pi, type = "h")



## close the raw file to clear memory
close(rawData)

## Shiny UI for peak and data summary
library(shiny)
# https://mastering-shiny.org/
#runExample("01_hello")
runApp("App-1")
