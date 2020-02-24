setwd("C:/Users/k1775475/OneDrive - King's College London/Work/Staff data/Mark data")

importData <- read.csv("GLA_LPOS_MSMS_features_17October2018.csv", check.names = FALSE)
metadata <- read.csv("GLA_simplified_database_17Oct2018.csv")

getOverlap <- intersect(importData$`GLA ID`, metadata$baseStudyNumber)
keepQCs <- importData[which(importData$`GLA ID` == "QC"),]

importData <- importData[which(is.element(importData$`GLA ID`, getOverlap)),]
metadata <- metadata[which(is.element(metadata$baseStudyNumber, getOverlap)),]
paste(importData$`GLA ID`) == paste(metadata$baseStudyNumber)

finalDataset <- cbind(importData[,c(1:8)], metadata, importData[,c(9:ncol(importData))])
importData <- finalDataset
importData <- importData[which(importDataset$XL == 0), ]
importData <- importData[which(importData$day == ""), ]
importData$patientGroup <- factor(importData$patientGroup)

## septic
importData$LDSseptic[which(importData$LDSseptic == "")] <- NA
tidySeptic <- as.numeric(importData$LDSseptic)
tidySeptic[which(tidySeptic == 3)] <- 1
tidySeptic[which(tidySeptic == 2)] <- 0
importData$LDSseptic <- tidySeptic

### subsetting and plotting
subset1 <- importData[which(importData$patientGroup == "ACLF"), ]
subset2 <- importData[which(importData$patientGroup == "AD"), ]
selectSubset <- rbind(subset1, subset2)
selectSubset$patientGroup <- factor(selectSubset$patientGroup)

removeCol <- which(colnames(selectSubset) == "fragment?")
removeCol <- which(colnames(selectSubset) == "Cyclopassifloic acid B")
selectSubset <- selectSubset[,-removeCol]


## produce plot

# survival 30 days
patientGroup <- selectSubset$survival30day
colNumbers <- c(404:ncol(selectSubset))
plotTitle <- "survival 30 day volcano plot putatively annotated features"

source("volcanoPlot.R")
getVolcano <- volcanoPlot(patientGroup, selectSubset, colNumbers, plotTitle) 

# survival 90 days
patientGroup <- selectSubset$survival90day
colNumbers <- c(404:ncol(selectSubset))
plotTitle <- "survival 90 day volcano plot putatively annotated features"

source("volcanoPlot.R")
getVolcano <- volcanoPlot(patientGroup, selectSubset, colNumbers, plotTitle) 

# survival 1 year
patientGroup <- selectSubset$survival1year
colNumbers <- c(404:ncol(selectSubset))
plotTitle <- "survival 1 year volcano plot putatively annotated features"

source("volcanoPlot.R")
getVolcano <- volcanoPlot(patientGroup, selectSubset, colNumbers, plotTitle) 

# sepsis
patientGroup <- selectSubset$LDSseptic
colNumbers <- c(404:ncol(selectSubset))
plotTitle <- "sepsis volcano plot putatively annotated features"

source("volcanoPlot.R")
getVolcano <- volcanoPlot(patientGroup, selectSubset, colNumbers, plotTitle) 
