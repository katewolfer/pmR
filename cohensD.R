##################################################
## Cohen's D function                           ##
## Kate Wolfer, Institute of Hepatology, London ##
## 06 December 2018                             ##
##################################################

## INPUT
## classifier: vector of classifier values, logical
## experimentalData: vector or dataframe of experimental data
##                   observations, in columns, numeric
## fileName: string for csv file output name

## OUTPUT
## dfClassifier: dataframe with feature names, mean and SD values
##               Cohen's D, Cohen's D significance as category 
##               (low, medium, high), Kruskal-Wallis test p-value
##               and Benjamini-Hochberg adjusted KW p-values
## csv file of dfClassifier
## only works for two-class case
## removes NA/NaN from the analysis to reduce potential dilution of classifier

cohensD <- function(classifier, experimentalData, fileName) {
  
  ## format classifier data
  classifier <- paste(classifier)
  getNA <- which(classifier == "NA")
  getNaN <- which(classifier == "NaN")
  classifier[which(classifier != "0")] <- "group 2"
  classifier[which(classifier == "0")] <- "group 1"
  classifier <- as.factor(classifier)
  
  # set up data frame variables
  SD <- NULL
  mean_N <- NULL
  mean_Y <- NULL
  KW_pval <- NULL
  Cohens_d <- NULL
  Cohens_sig <- NULL
  
  # loop through each experimental data column to collect required stats
  for (i in 1:ncol(experimentalData)){
    workingData <- data.frame(classifier, experimentalData[,i])
    workingData <- workingData[-c(getNA,getNaN),]
    statSD <- sd(workingData[,2])
    statMean <- tapply(workingData[,2],workingData$classifier, mean)
    SD[i] <- paste(statSD)
    mean_N[i] <- paste(statMean[1])
    mean_Y[i] <- paste(statMean[2])
    Cohens_d[i] <- (statMean[2] - statMean[1]) / statSD
    # assign "significance" for ease of interpretation
    if (Cohens_d[i] > 0.8 || Cohens_d[i] < -0.8){
      Cohens_sig[i] <- "high"
    } else if (Cohens_d[i] < 0.8 && Cohens_d[i] > 0.2 || Cohens_d[i] > -0.8 && Cohens_d[i] < -0.2){
      Cohens_sig[i] <- "medium"
    } else {Cohens_sig[i] <- "low"}
    # Kruskal-Wallis test - could be any other test, assumption of non-normality given that
    # we are looking at effect sizes between groups
    KWtest <- kruskal.test(experimentalData[,i] ~ classifier)
    KW_pval[i] <- KWtest$p.value
    rm(KWtest)
  }
  
  ## create result dataframe
  dfClassifier <- data.frame(paste(colnames(experimentalData)),mean_N, mean_Y, SD, Cohens_d, Cohens_sig, KW_pval)
  colnames( dfClassifier) <- c("feature","mean_0", "mean_1", "SD","Cohen's d","CD significance", "KW_pval")
  
  ## perform BH-adjustment of KW p-values
  dfClassifier$BH_KW <- p.adjust(dfClassifier$KW_pval, "BH")
  
  ## output csv file
  write.csv(dfClassifier, fileName)
  return(dfClassifier)
  
}



