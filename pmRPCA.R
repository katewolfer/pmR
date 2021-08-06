#' Title
#'
#' @param df
#' @param pc1
#' @param pc2
#' @param logTr
#' @param startCol
#' @param colorCol
#'
#' @return
#' @export
#'
#' @examples
#'
pmRPCA <- function(df, pc1, pc2, logTr, startCol, colorCol, labelPoints) {

  ## simple PCA, produces scores and loadings plots
  ## R package does not produce Q2 values based on 7-fold model
  ## cross-validation; we used SIMCA SIMCA+ 16.0 (Umetrics, Umea, Sweden)
  ## to produce the actual PCA, then used the ggplot2 code at the bottom
  ## of this function to produce the plots in the paper
  ## we will implement cross-validation at a later date

  ## INPUTS

  ## df: dataframe of collated atmospheric measurement data
  ## with columns as observations and rows as samples

  ## pc1, pc2: select which two principal components are plotted, this may need
  ## to be changed if PCs 1 and 2 do not show all of the interesting trends

  ## logTr: log10 transform data - mostly only suitable for LC-MS data

  ## startCol: integer of beginning of numeric observation data
  ## which needs to be normalised to ensure all values are actual numerical
  ## format, have standardised units no major outliers

  ## select column data used for coloring scores
  if (length(colorCol) == 0){
    colorData = 0
  } else if(length(colorCol) > 1){
    colorData <- colorCol
  } else if(length(colorCol) == 1){
    colorData <- df[,colorCol]
  }

  if (labelPoints == 0){
    labels <- 0
  } else if (length(labelPoints) == 1){
    labels <- df[,labelPoints]
  } else if (length(labelPoints) > 1){
    labels <- as.character(rownames(df))
  }



  ## select only the numerical data required to produce the model
  df <- df[,c(startCol:ncol(df))]

  # ## remove columns which are more than 50% missing values
  # pcmv <- NULL
  # for (cl in 1:ncol(df)){
  #   pcmv[cl] <- length(which(is.na(df[,cl])))/nrow(df)
  # }
  # forRemoval <- which(pcmv > 0.55)
  # df <- df[,-forRemoval]

  # ## unit variance scaling
  # getSD <- apply(df, 2, sd, na.rm = TRUE)
  # getMeans <- apply(df, 1, median, na.rm = TRUE)
  # for (cl in 1:ncol(df)){
  #   df[,cl] <- df[,cl]/getSD[cl]
  #   df[,cl] <- df[,cl]-getMeans[cl]
  # }

  ## PCA needs zero NA
  df[is.na(df)] <- 1e-19

  ## log10 transform data if required, produce PCA model
  set.seed(20)
  if (logTr == 1) {
    pcaSamples <- prcomp(log10(df+10))
    cat("data log10 transformed \n")
  } else if (logTr == 0) {
    pcaSamples <- prcomp(df)
    cat("data not transformed \n")
  } else (error("please specify if transformation to be applied: 0 for none, 1 for log(10)"))

  ## get PCA scores
  pcaDF <- as.data.frame(pcaSamples$x)
  pcaDF$group <- colorData

  ## get PCA loadings
  pcaLoadings <- as.data.frame(pcaSamples$rotation)
  # loadingsCol <- loadingsCol[c(startCol:ncol(df))]
  # loadingsCol <- loadingsCol[-forRemoval]
  # pcaLoadings$category <- loadingsCol
  # pcaLoadings$category <- as.factor(pcaLoadings$category)

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
  if (length(unique(colorData)) < 5){
  } else
    if(is.numeric(colorData)){
      scoresPlot <- scoresPlot + scale_color_viridis_c(option = "plasma")
    } else {
      scoresPlot <- scoresPlot + scale_color_viridis_d(option = "plasma")
    }
  # scoresPlot <- scoresPlot + ggtitle(paste0("PCA scores plot, points coloured by ",
  #                                           colnames(df)[colorCol]))

  if (labels != 0){
  scoresPlot <- scoresPlot + geom_label_repel(aes(label = labels))
  }

  # forceColors <- rainbow(length(unique(labels)))
  # forceColors <- c("red", "red","blue","blue", "orange","orange",
  #                  "yellow", "yellow", "pink", "pink",
  #                  "grey","grey","grey","grey",
  #                  "grey","grey","grey","grey",
  #                  "black")
  #forceColors <- rainbow(length(unique(labels)))
  #forceColors <- c("#241984", "magenta","#ff8600","#fff105", "#6633ff")
  #scoresPlot <- scoresPlot + scale_color_manual(values = forceColors)

  scoresPlot <- scoresPlot + geom_point(data = getHT2, aes(x = V1, y = V2),
                                        colour = '#888888', size = 0.2)
  scoresPlot <- scoresPlot + geom_hline(yintercept = 0)
  scoresPlot <- scoresPlot + geom_vline(xintercept = 0) + theme_bw()
  scoresPlot <- scoresPlot + theme(text = element_text(size = 16))
  scoresPlot <- scoresPlot + xlab(paste0("PC", pc1, "\n R2X ", expl.var[pc1], "%"))
  scoresPlot <- scoresPlot + ylab(paste0("PC", pc2, "\n R2X ", expl.var[pc2], "%"))
  # scoresPlot <- scoresPlot + theme(plot.title = element_text(hjust = 0.5),
  #                                  legend.title=element_blank())


  ## sort the column names
  loadLabel <- gsub("X", "", colnames(df))

  ## obtain loadings plot
  loadingsPlot <- ggplot(pcaLoadings, aes(x = pcaLoadings[,pc1],
                                          y = pcaLoadings[,pc2]))
  loadingsPlot <- loadingsPlot + geom_hline(yintercept = 0)
  loadingsPlot <- loadingsPlot + geom_vline(xintercept = 0) + theme_bw()
  loadingsPlot <- loadingsPlot + geom_point(size = 5.5, colour = "blue")
  loadingsPlot <- loadingsPlot + ggtitle(paste0("PCA loadings plot"))
  # loadingsPlot <- loadingsPlot + geom_text(aes(label=loadLabel),
  #                                          hjust=-0.25, color = "#888888",
  #                                          size = 3.5)


  loadingsPlot <- loadingsPlot + xlab(paste0("PC", pc1, "\n R2X ",
                                             expl.var[pc1], "%"))
  loadingsPlot <- loadingsPlot + ylab(paste0("PC", pc2, "\n R2X ",
                                             expl.var[pc2], "%"))
  loadingsPlot <- loadingsPlot + theme(text = element_text(size = 16))
  # loadingsPlot <- loadingsPlot + theme(legend.title=element_blank())
  # loadingsPlot <- loadingsPlot + theme(legend.position="bottom")
  loadingsPlot <- loadingsPlot + theme(plot.title = element_text(hjust = 0.5))
  loadingsPlot

  return(list(scoresPlot, loadingsPlot, pcaDF))

}
