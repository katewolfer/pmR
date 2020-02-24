################################
## Visualising scaling in MVA ##
## Kate Leary                 ##
## 08 January 2020            ##
################################

# load up libraries
library(ggplot2)

# import data and rearrange
dataImport <- read.csv("test_data.csv") # import data, change file title as needed
msData <- as.data.frame(t(dataImport[,11:ncol(dataImport)])) # change 11 if needed to reflect number of non-feature columns
colnames(msData) <- dataImport[,1]

## unscaled data - collect means and SD for each column
collectMeans <- apply(msData[,1:ncol(msData)],2,mean)
collectSDs <- apply(msData[,1:ncol(msData)],2,sd)

# sort the SDs by mean intensity rank
getRank <- sort(collectMeans, index.return=TRUE)
rankSDs <- collectSDs[getRank$ix]

# collect the rank number and SD value into a data frame for plotting
plotting <- as.data.frame(cbind(1:ncol(msData),rankSDs))
colnames(plotting)[1] <- "featureRank"


## Plotting, can be adjusted
png("rank intensity SD plot.png", h=2500, w=4000, res=400) #change title if needed for new plot
p1 <- ggplot(plotting, aes(x=plotting$featureRank, y=plotting$rankSDs)) + 
  geom_point(alpha=0.75,color = '#3333FF') 
p1 <- p1 + theme_bw()
#p1 <- p1 + theme(legend.position="none")
#p1 <- p1 + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 16)) # if you need a title - put it in between the ""
#p1 <- p1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) # if you want to remove the spare space
p1 <- p1 + xlab("rank of mean intensity") 
p1 <- p1 + ylab("standard deviation")
p1
dev.off()

#######################
## Comparing scaling ##
#######################

# pareto scaled data
# 1/(Sk)1/2, where Sk is the standard deviation
paretoData <- msData
for (i in 1:ncol(msData)) {
  paretoData[,i] <- msData[,i] /(1/sqrt(sd(msData[,i]))) }

paretoMean <- msData
for (i in 1:ncol(msData)) {
  paretoMean[,i] <- (msData[,i]/(mean(msData[,i])))/(1/sqrt(sd(msData[,i]))) }

# mean-centered data
meanCentered <- msData
for (i in 1:ncol(msData)) {
  meanCentered[,i] <- msData[,i]/(mean(msData[,i])) }

# setting up plotting
matchColours = rainbow(ncol(dataImport))

maxMS <- max(msData[,1:5])
maxMean <- max(meanCentered[,1:5])
maxPar <- max(paretoData[,1:5])
minMS <- min(msData[,1:5])
minMean <- min(meanCentered[,1:5])
minPar <- min(paretoData[,1:5])

# make plot
png("scaling plot.png", h=5000, w=3000, res=300) #change title if needed for new plot
par(mfrow=c(5,3))
for (i in 1:5) {
  plot(msData[,i], ylab = "untransformed data", col = matchColours, pch = 19, ylim = c(minMS,maxMS))
  plot(meanCentered[,i], ylab = "mean-centered data",col = matchColours, pch = 19, ylim = c(minMean,maxMean))
  plot(paretoData[,i], ylab = "pareto-scaled data",col = matchColours, pch = 19, ylim = c(minPar,maxPar)) }
dev.off()

png("log scaling plot.png", h=5000, w=3000, res=300) #change title if needed for new plot
par(mfrow=c(5,3))
for (i in 1:5) {
  plot(log10(msData[,i]+10), ylab = "untransformed data", col = matchColours, pch = 19, ylim = c(log10(minMS+10),log10(maxMS+10)))
  plot(log10(meanCentered[,i]+10), ylab = "mean-centered data",col = matchColours, pch = 19, ylim = c(log10(minMean+10),log10(maxMean+10)))
  plot(log10(paretoData[,i]+10), ylab = "pareto-scaled data",col = matchColours, pch = 19, ylim = c(log10(minPar+10),log10(maxPar+10))) }
dev.off()

## second comparison
maxMean <- max(meanCentered[,1:5])
maxParM <- max(paretoMean[,1:5])
maxPar <- max(paretoData[,1:5])
minMean <- min(meanCentered[,1:5])
minParM <- min(paretoMean[,1:5])
minPar <- min(paretoData[,1:5])

# make plot
png("scaling plot 2.png", h=5000, w=3000, res=300) #change title if needed for new plot
par(mfrow=c(5,3))
for (i in 1:5) {
  plot(meanCentered[,i], ylab = "mean-centered data",col = matchColours, pch = 19, ylim = c(minMean,maxMean))
  plot(paretoMean[,i], ylab = "pareto-scaled centered data",col = matchColours, pch = 19, ylim = c(minParM,maxParM)) 
  plot(paretoData[,i], ylab = "pareto-scaled data",col = matchColours, pch = 19, ylim = c(minPar,maxPar)) }
dev.off()

png("log scaling plot 2.png", h=5000, w=3000, res=300) #change title if needed for new plot
par(mfrow=c(5,3))
for (i in 1:5) {
  plot(log10(meanCentered[,i]+10), ylab = "mean-centered data",col = matchColours, pch = 19, ylim = c(log10(minMean+10),log10(maxMean+10)))
  plot(log10(paretoMean[,i]+10), ylab = "pareto-scaled centered data",col = matchColours, pch = 19, ylim = c(log10(minParM+10),log10(maxParM+10))) 
  plot(log10(paretoData[,i]+10), ylab = "pareto-scaled data",col = matchColours, pch = 19, ylim = c(log10(minPar+10),log10(maxPar+10))) }
dev.off()


