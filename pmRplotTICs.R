pmRplotTICs <- function(){

  require(MSnbase)
  require(ggplot2)

  ## list all the mzML files in the directory
  fileList <- list.files(getwd(), pattern = "mzML")

  ## pull out all the BPI data
  raw_data <- readMSData(files = fileList,mode = "onDisk")
  chr_raw <- chromatogram(raw_data,aggregationFun = "max")
  cat("All raw data imported \n")

  ## reshape the data as Xcalibur saves the data with time disjoints
  pullRT <- NULL
  for(b in 1:length(chr_raw)){
    pullRT <- c(pullRT, chr_raw@.Data[[b]]@rtime)
  }

  pullRT <- unique(sort(pullRT))
  df <- as.data.frame(pullRT)
  rownames(df) <- c(1:nrow(df))
  colnames(df) <- "rt"

  for(k in 1:length(chr_raw)){
    for(c in 1:nrow(df)){
      findRT <- which(df[c,1] == chr_raw@.Data[[k]]@rtime)
      if(length(findRT) > 0){
        df[c,k+1] <- chr_raw@.Data[[k]]@intensity[findRT]
      }
    }
  }
  cat("All RT points aligned with intensity values \n")
  colnames(df)[c(2:(length(chr_raw)+1))] <- chr_raw@phenoData[[1]]

  ## melt the dataframe because ggplot2 prefers this
  meltFeature <- melt(df, id.vars = "rt")
  ## this line needed to get actual line plot, else loess smoothing needed
  ## loess span does not handle this level of resolution well
  meltFeature <- meltFeature[-which(is.na(meltFeature$value)),]

  ## construct the plot
  TICplot <- ggplot(data = meltFeature, aes(x = rt, y = value, color = variable))
  TICplot <- TICplot + geom_line(size = 0.7, na.rm = TRUE)
  TICplot <- TICplot + xlab('retention time (mins)') + ylab('intensity (arbitrary units)')
  TICplot <- TICplot + theme(plot.title = element_text(hjust = 0.5, size = 18))
  TICplot <- TICplot + theme(axis.line = element_line(colour = "black"),
                             panel.grid.minor = element_line(colour="light grey",
                                                             size=0.01),
                             panel.border = element_blank(),
                             panel.grid.major = element_line(colour="light grey",
                                                             size=0.01),
                             panel.background = element_blank())
  TICplot <- TICplot + scale_x_continuous(expand = c(0,0))
  TICplot <- TICplot + scale_y_continuous(expand = c(0,0))
  TICplot <- TICplot + labs(colour='sample') + theme(legend.position = c(0.90,0.90))

  return(TICplot)

}


