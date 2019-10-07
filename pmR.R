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
require(rgl)
require(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("xcms", lib="C:/Program Files/R/R-3.5.3/library")
BiocManager::install("BiocParallel", lib="C:/Program Files/R/R-3.5.3/library")


## set Proteowizard program path
proteowizardPath <- 'C:/Program Files/ProteoWizard/ProteoWizard 3.0.19276.7124c5404'

## set target folder where raw files live
targetFolder     <- normalizePath(file.path('C:', 
                                            'Users', 
                                            'Kate', 
                                            'Desktop',
                                            'OrbitrapData'))

# create new directory for mzML converted data
dir.create(paste0(targetFolder, '\\converted_data'))
saveFolder        <- paste0(targetFolder, '\\converted_data')

# perform file conversion
source("pmRfileConvert.R")
pmRfileConvert(proteowizardPath, targetFolder, saveFolder)

# set the working directory to the location of the converted files
setwd(saveFolder)

# parameters
params <- pmRparameterization(getwd(), 10) 


# top50 <- sort(checkFile@tic, decreasing = TRUE)
# summary(top50[1:(round(length(checkFile@tic)/2))])
# summary(top50[round(length(checkFile@tic)/2):length(top50)])

checkFile1 <- xcmsRaw(fileList[1], profstep = 0.1)
checkFile2 <- xcmsRaw(fileList[2], profstep = 0.1)
checkFile3 <- xcmsRaw(fileList[3], profstep = 0.1)
checkFile4 <- xcmsRaw(fileList[4], profstep = 0.1)

## determine experiment polarity from file data
## need to add check for polarity switching files
checkNeg <- sum(checkFile@polarity=='negative')
checkPos <- sum(checkFile@polarity=='positive')
if (checkNeg > checkPos){
  polarity = 'negative'
} else { polarity ='positive'}

## how to make a BPI chromatogram?
plot(checkFile@scantime/60, 
     checkFile@tic, 
     type = 'l', 
     col = 'blue',
     xlab = 'RT (mins)',
     ylab = paste0('intensity, ', polarity, ' ESI mode'))

## plot all TICs
allTICs = cbind(checkFile1@scantime[c(1:692)]/60,
                checkFile1@tic[c(1:692)], 
                checkFile2@tic[c(1:692)], 
                checkFile3@tic[c(1:692)], 
                checkFile4@tic[c(1:692)])
allTICs = as.data.frame(allTICs)
colnames(allTICs) = c('RT', 'TIC_1', 'TIC_2','TIC_3','TIC_4')

plot(allTICs$RT, 
     allTICs$TIC_1,
     type = 'l', 
     col = 'blue',
     xlab = 'RT (mins)',
     ylab = paste0('intensity, ', polarity, ' ESI mode'))

par(new=TRUE)
plot(allTICs$RT, 
     allTICs$TIC_2,
     type = 'l', 
     col = 'dark green',     
     xlab = '',
     ylab = '')

par(new=TRUE)
plot(allTICs$RT, 
     allTICs$TIC_3,
     type = 'l', 
     col = 'red',     
     xlab = '',
     ylab = '')

par(new=TRUE)
plot(allTICs$RT, 
     allTICs$TIC_4,
     type = 'l', 
     col = 'purple',     
     xlab = '',
     ylab = '')


getCols <- brewer.pal(n = 10, name = "Dark2")

p1 <- ggplot(data = allTICs)
for (i in c(2:5)) {
  p1 <- p1 + geom_line(data = allTICs, 
                       aes(x=allTICs[,1],
                           y=allTICs[,i]), 
                       colour = getCols[i], 
                       size = 0.25) }
p1 <- p1 +  xlab('RT (mins)') + ylab('intensity')
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5, size = 18))
p1 <- p1 + theme(axis.line = element_line(colour = "black"),
                 panel.grid.minor = element_line(colour="light grey", size=0.01),
                 panel.border = element_blank(),
                 panel.grid.major = element_line(colour="light grey", size=0.01),
                 panel.background = element_blank())
p1 <- p1 + scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(allTICs[,1]), by = 5))
p1 <- p1 + scale_y_continuous(expand = c(0, 0))
p1



cdffiles <- list.files(getwd(), recursive = TRUE, full.names = TRUE)
getTICs(files=cdffiles, pdfname="fileTICs.pdf")
