pmRQCbatchChecking <- function(){

  require(plotly)
  require(xcms)

  ## data extraction
  fileList <- list.files(getwd(), pattern = "mzML")
  xs <- xcmsSet(fileList[grep("QC", fileList)],
                method = "centWave",
                ppm = 2.5,
                BPPARAM = SnowParam(4))

  bw = 6
  MZerror = 0.005
  gagds <- group(xs, method="density", mzwid=MZerror, bw=bw)
  fds <- fillPeaks(gagds)
  fill <- fds
  vals.raw <- groupval(fill, "medret", "into")
  vals.raw.meta <- cbind(groups(fds), vals.raw)
  df <- as.data.frame(t(vals.raw))

  pc1 <- 1
  pc2 <- 2

  ## PCA needs zero NA
  df[is.na(df)] <- 1e-19

  ## log10 transform data if required, produce PCA model
  set.seed(20)
  #if (logTr == 1) {
    pcaSamples <- prcomp(log10(df+10))
    cat("data log10 transformed \n")
  #} else if (logTr == 0) {
  #  pcaSamples <- prcomp(df)
  #  cat("data not transformed \n")
  #} else (error("please specify if transformation to be applied: 0 for none, 1 for log(10)"))

  ## get PCA scores
  pcaDF <- as.data.frame(pcaSamples$x)
  pcaDF$group <- c(1:nrow(pcaDF))

  ## get PCA loadings
  pcaLoadings <- as.data.frame(pcaSamples$rotation)

  ## R2X value
  expl.var <- round(pcaSamples$sdev^2/sum(pcaSamples$sdev^2)*100,2)

  ## get Hotellings T2 value for plot
  getHT2 <- as.data.frame(pcaMethods:::simpleEllipse(pcaDF[,"PC1"],
                                                     pcaDF[,"PC2"],
                                                     alfa = 0.95,
                                                     len = 500))

  ## obtain plot
  scoresPlot <- ggplot(pcaDF, aes(x = pcaDF[,pc1], y = pcaDF[,pc2],
                                  color = group)) +
    geom_point(size = 5, alpha = 0.8)

  scoresPlot <- scoresPlot + scale_color_viridis_c(option = "plasma")

  scoresPlot <- scoresPlot + geom_point(data = getHT2, aes(x = V1, y = V2),
                                        colour = '#888888', size = 0.2)
  scoresPlot <- scoresPlot + geom_hline(yintercept = 0)
  scoresPlot <- scoresPlot + geom_vline(xintercept = 0) + theme_bw()
  scoresPlot <- scoresPlot + theme(text = element_text(size = 16))
  scoresPlot <- scoresPlot + xlab(paste0("PC", pc1, "\n R2X ", expl.var[pc1], "%"))
  scoresPlot <- scoresPlot + ylab(paste0("PC", pc2, "\n R2X ", expl.var[pc2], "%"))


  ## sort the column names
  loadLabel <- gsub("X", "", colnames(df))

  ## obtain loadings plot
  loadingsPlot <- ggplot(pcaLoadings, aes(x = pcaLoadings[,pc1],
                                          y = pcaLoadings[,pc2]))
  loadingsPlot <- loadingsPlot + geom_hline(yintercept = 0)
  loadingsPlot <- loadingsPlot + geom_vline(xintercept = 0) + theme_bw()
  loadingsPlot <- loadingsPlot + geom_point(size = 5.5, colour = "blue")
  loadingsPlot <- loadingsPlot + ggtitle(paste0("PCA loadings plot"))


  loadingsPlot <- loadingsPlot + xlab(paste0("PC", pc1, "\n R2X ",
                                             expl.var[pc1], "%"))
  loadingsPlot <- loadingsPlot + ylab(paste0("PC", pc2, "\n R2X ",
                                             expl.var[pc2], "%"))
  loadingsPlot <- loadingsPlot + theme(text = element_text(size = 16))
  loadingsPlot <- loadingsPlot + theme(plot.title = element_text(hjust = 0.5))
  loadingsPlot

  ## make a 3D PCA plot
  stitchData <- as.data.frame(pcaDF)
  #stitchData$sample <- df[,5]
  stitchData$sample <- c(1:nrow(stitchData))


  p <- plot_ly(stitchData,
               x = stitchData$PC1,
               y = stitchData$PC2,
               z = stitchData$PC3,
               color = stitchData$sample,
               colors = c("#241984", "magenta","#fff105")
  ) %>%
    add_markers() %>%
    layout(
      scene = list(xaxis = list(title = 'PC1'),
                   yaxis = list(title = 'PC2'),
                   zaxis = list(title = 'PC3'))
    )

  p

  return(list(pcaDF, loadingsPlot, scoresPlot, p))
}





