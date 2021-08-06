pmRmanualFeatures <- function(filterTable, fileList, raw_data){

  require(MSnbase)
  require(ggplot2)
  require(reshape2)

  # raw_data <- readMSData(files = fileList,
  #                        mode = "onDisk")

  ## add columns to the table for annotation
  filterTable$realFeature <- NA
  filterTable$notes <- NA

  ## systematically run through each feature, plot, then
  ## ask user for manual input on retaining peak in final
  ## data table, export figure to image file if kept
  for(i in 1:nrow(filterTable)){

    ## output for progress
    cat(paste0("feature ", i, "/",nrow(filterTable)))
    cat("\n")

    # MSnbase extract feature params
    rtr <- c(filterTable$rtmed[i]-20,filterTable$rtmed[i]+20)
    mzr <- c(filterTable$mzmed[i]-0.0006, filterTable$mzmed[i]+0.0006)

    ## extract the chromatogram
    chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)

    # ## collect useful data about the feature
    # peakData <- data.frame(min = NA,
    #                        max = NA,
    #                        skew = NA,
    #                        kurtosis = NA,
    #                        se = NA,
    #                        stringsAsFactors=FALSE)

    # ## assess the peaks
    # for (c in 1:length(chr_raw)){
    #   ## sometimes the data is too jaggy to plot line so use loess
    #   getPeakSmooth <- loess.smooth(x = chr_raw@.Data[[c]]@rtime,
    #                                 y =  chr_raw@.Data[[c]]@intensity,
    #                                 span = 0.15, degree = 1,
    #                                 family = "gaussian",
    #                                 xlab = NULL, ylab = NULL)
    #   plot(getPeakSmooth)
    #   readline("?")
    #   peakDescriptors <- describe(getPeakSmooth$y)
    #   peakData[c,] <- unlist(peakDescriptors[c(8,9,11,12,13)])
    # }

    ## make ggplot2, plot raw and smoothed side by side?
    pullRT <- NULL
    for(b in 1:length(chr_raw)){
      pullRT <- c(pullRT, chr_raw@.Data[[b]]@rtime)
    }

    pullRT <- unique(sort(pullRT))

    df <- as.data.frame(pullRT)
    rownames(df) <- c(1:nrow(df))
    colnames(df) <- "rt"

    getTitle <- paste0("Feature ", i, ", ", round(filterTable$rtmed[i]/60,2),
                       " mins, mz ",round(filterTable$mzmed[i],5))
    cat(getTitle)
    cat("\n")

    for(k in 1:length(chr_raw)){
      for(c in 1:nrow(df)){
        for(j in 1:length(chr_raw@.Data[[k]]@intensity))
          findRT <- which(df[c,1] == chr_raw@.Data[[k]]@rtime)
        #cat("found RT")
        if(length(findRT) > 0){
          df[c,k+1] <- chr_raw@.Data[[k]]@intensity[findRT]
          #cat("found m/z")
        } #else {cat("m/z not found!")}
        #cat("done!")
      }

    }

    colnames(df)[c(2:(length(chr_raw)+1))] <- chr_raw@phenoData[[1]]

    meltFeature <- melt(df, id.vars = "rt")

    p5 <- ggplot(data = meltFeature, aes(x = rt, y = value, color = variable))
    p5 <- p5 + geom_point(size = 0.7, na.rm = TRUE)
    p5 <- p5 + geom_smooth(method="loess", se=F, span = 0.1, size = 0.2)
    p5 <- p5 + xlab('retention time (mins)') + ylab('intensity (arbitrary units)')
    p5 <- p5 + theme(plot.title = element_text(hjust = 0.5, size = 18))
    p5 <- p5 + theme(axis.line = element_line(colour = "black"),
                     panel.grid.minor = element_line(colour="light grey", size=0.01),
                     panel.border = element_blank(),
                     panel.grid.major = element_line(colour="light grey", size=0.01),
                     panel.background = element_blank())
    p5 <- p5 + scale_x_continuous(expand = c(0,0))
    p5 <- p5 + scale_y_continuous(expand = c(0,0))
    p5 <- p5 + labs(colour='sample') + theme(legend.position = c(0.85,0.85))

    p5 <- p5 + ggtitle(getTitle)
    plot(p5)


    ## annotate
    checkFeature <- readline(prompt="Seems like a real feature?: ")
    checkFeature <- tolower(checkFeature)

    ## if Y, save the plot in publication quality, open as PDF
    if (checkFeature == "y"){
      filterTable$realFeature[i] <- "Y"
      ggsave(paste0(getTitle,".eps"),
             plot = p5,
             device = "eps",
             path = NULL,
             scale = 1,
             width = 12.5,
             height = 8,
             dpi = 600)
      
      ggsave(paste0(getTitle,".png"),
             plot = p5,
             device = "png",
             path = NULL,
             scale = 1,
             width = 12.5,
             height = 8,
             dpi = 600)

             ## automated labelling for time saving
             cat("Choose from the following annotations: \n")
             cat("1: select for further investigation \n")
             cat("2: duplicated feature, still of interest \n")
             cat("3: low intensity but retain for checking \n")
             getClass <- readline(prompt="Enter label category: ")
             if(getClass == 1){
               getNotes <- "further investigation"
             } else if(getClass == 2){
               getNotes <- "duplicated feature Y"
             } else if(getClass == 3){
               getNotes <- "low intensity"
             } else {getNotes <- readline(prompt="Notes about the peak: ")}

    } else {filterTable$realFeature[i] <- "N"
    cat("Choose from the following annotations; any other category can have manual notes added: \n")
    cat("1: looks like baseline/poorly resolved \n")
    cat("2: duplicated feature, not of interest \n")
    cat("3: too low intensity \n")
    getClass <- readline(prompt="Enter label category: ")
    if(getClass == 1){
      getNotes <- "baseline/poorly resolved N"
    } else if(getClass == 2){
      getNotes <- "duplicated feature"
    } else if(getClass == 3){
      getNotes <- "too low intensity"
    } else {getNotes <- readline(prompt="Notes about the peak: ")}
    }

    filterTable$notes[i] <- getNotes
    dev.off()
  }

  return(filterTable)
}


