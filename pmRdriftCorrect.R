pmRdriftCorrect <- function(ds, startCol, batchIndex, curveSpan){

  toRemove <- c(which(ds$runType == "blank"),
                which(ds$runType == "dilution"),
                which(ds$runType == "conditioning"))
  SRs <- which(ds$runType[-toRemove] == "QC")
  x <- as.numeric(paste(ds$absRunOrder[-toRemove]))
  df <- ds[-toRemove,c(startCol:ncol(ds))]

  ## select 50 random features for correction plotting
  sampleRand <- sort(sample(1:ncol(df),50))

  ## set up container for corrected data
  corrected.ds <- NULL

  ## find non-zero minimum
  collectMin <- NULL
  for(z in 1:ncol(df)){
    collectMin[z] <- min(df[(df[,z] > 0),z])
  }
  fixMin <- min(collectMin)

  ## iterate through features and apply correction
  for (i in 1:ncol(df)){

    ## select feature
    target <- i
    print(paste0("Drift correcting feature #", target))

    ## blank dummy vector
    y <- rep(NA, length(x))

    ## fills in values ONLY for SR sample
    y[SRs] <- as.numeric(df[SRs,target])

    ## alternative NA fill
    # if(is.na(y[length(y)])){
    #   boundY <- length(y+1)
    #   y[length(y)+1] <- as.numeric(df[max(SRs),target])
    #   x[length(x)+1] <- max(x)+1
    # } else if (is.na(y[1])){
    #   boundY <- 1
    #   y[1] <- as.numeric(df[min(SRs),target])
    #   x <- c(min(x)-1, x)
    # }

    ## add upper and lower boundaries for absolute permitted corretion values
    upper.lim <- median(y, na.rm=T)+(2*sd(y, na.rm=T))
    lower.lim <- median(y, na.rm=T)-(2*sd(y, na.rm=T))
    y[which(y > upper.lim | y < lower.lim)] <- NA


    # fits the loess curve to the SR data only
    y.loess <- stats::loess(y ~ x,
                            span=spam,
                            data.frame(x=x, y=y),
                            control=loess.control(surface='interpolate'))

    # generates predicted values for all x positions (all samples)
    y.predict <- predict(y.loess, data.frame(x=x))

    # caps predicted vals at the min/max in the original data
    y.predict[which(y.predict>max(y[SRs]))]<- max(y[SRs])
    y.predict[which(y.predict<min(y[SRs]))]<- min(y[SRs])

    ## if correction has not run from beginning/until end of the vector
    bracket.check.end <- is.na(y.predict[nrow(df)])
    if (bracket.check.end) {
      # finds last non NA value (sets as k)
      for (k in nrow(df):(0.5*nrow(df))){
        end.check <- is.na(y.predict[k])
        if (end.check) next else break
      }
      y.predict[k:nrow(df)] <- y.predict[k]
    } else NULL

    bracket.check.start <- is.na(y.predict[1])
    if (bracket.check.start) {
      for (k in 1:(0.5*nrow(df))){
        start.check <- is.na(y.predict[k])
        if (start.check) next else break
      }
      y.predict[1:k] <- y.predict[k]
    } else NULL

    ## summary of feature
    getSummary <- paste(summary(df[,i]))

    ## remove added dummy values if alternative NA fill used
    # if(boundY == 1){
    #   y <- y[-1]
    #   x <- x[-1]
    #   y.predict <- y.predict[-1]
    # } else {
    #   y <- y[-length(y)]
    #   x <- x[-length(x)]
    #   y.predict <- y.predict[-length(y.predict)]
    # }

    ## generates a correction vector that specifies multiplication factor
    ## for each sample to be adjusted up to max predicted value
    correction <- max(y.predict, na.rm=T)/y.predict
    correction.applied <- as.numeric(df[,target])*correction
    correction.applied[which(correction.applied < as.numeric(getSummary[1]))] <- as.numeric(getSummary[1])

    ## binds the corrected vector to the growing matrix of corrected feature values
    corrected.ds <- rbind(corrected.ds, correction.applied)
    corrected.ds[i,] <- correction.applied

    ## for preselected random samples - plot before and after the correction
    if (is.element(i, sampleRand)) {
      ## stitch the data together for plotting
      p <- data.frame(x,
                      as.numeric(df[,target]),
                      y.predict,
                      correction.applied,
                      ds$runType[-toRemove],#[selectBatch],
                      ds$absRunOrder[-toRemove])#[selectBatch])

      # length(x)
      # length(as.numeric(df[,target]))
      # length(y.predict)
      # length(correction.applied)
      # length(ds$runType[-toRemove])
      # length(ds$absRunOrder[-toRemove])
      colnames(p) <- c("x", "it", "predict", "corr.it", "type","run.order")

      ## get maximum value to ensure plots have identical Y-axis
      maxY <- (max(correction.applied))*1.01
      #maxY <- 1e8

      ## before correction plot
      dp1 <- ggplot(p, aes(x=x, y=it, col=type)) +
        geom_point(shape=19, size=3,alpha=0.6) +
        geom_line(aes(x=x, y=predict)) +
        geom_text(data = p, aes(x = x, y = it, label = run.order),
                  hjust = 0, nudge_x = 0.5) +
        scale_colour_brewer(palette = "Set1") +
        labs(x="run order", y="feature intensity") +
        ggtitle(paste("feature", i, "(pre drift correction)", sep=" ")) +
        theme_bw()
      dp1 <- dp1 + coord_cartesian(ylim=c(0,maxY))
      dp1 <- dp1 + scale_y_continuous(expand = c(0, 0))

      ## after corretion plot
      dp2 <- ggplot(p, aes(x=x, y=corr.it, col=type)) +
        geom_point(shape=19, size=3,alpha=0.6) +
        geom_text(data = p, aes(x = x, y = corr.it, label = run.order),
                  hjust = 0, nudge_x = 0.5) +
        scale_colour_brewer(palette = "Set1") +
        labs(x="run order", y="feature intensity") +
        ggtitle(paste("feature", i, "(post drift correction)", sep=" ")) +
        theme_bw()
      dp2 <- dp2 + coord_cartesian(ylim=c(0,maxY))
      dp2 <- dp2 + scale_y_continuous(expand = c(0, 0))

      ## output plots for manual checking
      getPlot <- ggarrange(dp1, dp2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

      ## save plot to PNG file
      ggsave(file=paste(i," correction, batch", batchIndex, " ", Sys.Date(),".png", sep=""), plot=getPlot, height=15, width=45, units="cm", dpi=300)

    }
  }

  rownames(corrected.ds) <- NULL
  return(corrected.ds)
}
