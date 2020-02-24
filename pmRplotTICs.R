pmRplotTICs <- function(fileList){

  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################

  ## function to automatically plot TICs in the directory of converted files

  if(exists(fileList)){
  } else{
    fileList <- list.files(getwd())
  }

  getPolarity <- NULL
  lengthTIC <- NULL
  lengthSC <- NULL

  for (i in 1:length(fileList)){
    checkFile <- xcmsRaw(fileList[i], profstep = 10)
    if(i == 1){
      allTICs <- list(checkFile@tic)
      allSC <- list(checkFile@scantime)

    } else {
      allTICs[i] <- list(checkFile@tic)
      allSC[i] <- list(checkFile@scantime)
    }

    lengthTIC[i] <- length(checkFile@tic)
    lengthSC[i] <- length(checkFile@scantime)
    getPolarity[i] <- levels(factor(checkFile@polarity))
  }

  ## define the run polarity - this is not a failsafe method for other people's
  ## sample lists!!!
  checkPos <- sum(getPolarity == 'positive')
  checkNeg <- sum(getPolarity == 'negative')

  if (checkPos == 0){
    definePolarity = 'positive'
  } else if (checkNeg == 0){
    definePolarity ='positive'
  } else if (checkNeg > checkPos){
    definePolarity = 'negative'
  } else if (checkPos > checkNeg){
    definePolarity ='positive'}

  cat(paste0("Polarity of run: ", definePolarity))

  ## get minimum scan length for all TICs
  alignTimes <- NULL
  for (i in 1:length(allSC)){
    alignTimes[i] <- length(unlist(allSC[i]))
  }
  minSC <- min(alignTimes)

  ## scantimes
  checkSC <- unlist(allSC[1])
  checkSC <- as.matrix(checkSC[1:minSC])
  for(i in 2:length(allSC)){
    pullSC <- as.matrix(unlist(allSC[i]))
    checkSC <- cbind(checkSC, pullSC[1:minSC])
  }

  cat("Collecting scantimes ...")

  ## retention times
  checkTICs <- unlist(allTICs[1])
  checkTICs <- as.matrix(checkTICs[1:minSC])
  for(i in 2:length(allTICs)){
    pullTICs <- as.matrix(unlist(allTICs[i]))
    checkTICs <- cbind(checkTICs, pullTICs[1:minSC])
  }

  cat("Collecting TIC values ...")

  allTICs = as.data.frame(cbind(checkSC[,1],checkTICs))
  colnames(allTICs) <- c('RT', fileList)

  ## melt dataframe for faster plotting
  meltTICs <- melt(allTICs, id.vars = "RT")
  findMax <- max(meltTICs$value, na.rm = TRUE)
  meltTICs$RT <- as.numeric(paste(meltTICs$RT))/60

  ## produce single plot of all TICs
  p1 <- ggplot(meltTICs, aes(x=RT, y=value, col=variable)) + geom_line()
  p1 <- p1 +  xlab('RT (mins)') + ylab('intensity')
  p1 <- p1 + theme(plot.title = element_text(hjust = 0.5, size = 18))
  p1 <- p1 + theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour="light grey", size=0.01),
                   panel.border = element_blank(),
                   panel.grid.major = element_line(colour="light grey", size=0.01),
                   panel.background = element_blank())
  p1 <- p1 + scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(allTICs[,1]), by = 1))
  p1 <- p1 + scale_y_continuous(expand = c(0, 0), breaks = seq(0, findMax, by = 20000000))
  p1 <- p1 + labs(color=paste(definePolarity, 'ESI mode chromatograms'))

  return(p1)

}
