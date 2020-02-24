##################################################
## Tidy chromatogram plotting function          ##
## Kate Wolfer, Institute of Hepatology, London ##
## 06 December 2018                             ##
##################################################

## INPUT
## chrom: csv file, with two columns: RTmins, showing time, and
## intensity, showing intensity values for each timepoint
## cutStart: time in mins before which data should be discarded
## cutEnd: time in mins after which data should be discarded
## zoom: logical, 1 if subsection is required
## zoomSection: vector (e.g c(1,5)) indicating start and end of 
## RTmins to output zoom section of chromatogram

## OUTPUT

plotChromatogram <- function(chrom, cutStart, cutEnd, fileName, zoom, zoomSection) {
  
  require(ggplot2)
  require(RColorBrewer)
  
  ## remove any unwanted data e.g. cut first two mins and last five mins
  chrom <- chrom[which(chrom$RTmins < cutStart),]
  chrom <- chrom[which(chrom$RTmins > cutEnd),]
  
  ## if layered chromatograms are required
  getCols <- brewer.pal(n = 10, name = "RdBu")
  
  ## create plot
  p1 <- ggplot(data = chrom, aes(x=chrom$RTmins,y=chrom$sample1))
  p1 <- p1 + geom_line(size = 0.25, colour = getCols[1])
  
  if (ncol(chrom) > 2) {
    for (i in 3:ncol(chrom)) {
      p1 <- p1 + geom_line(data = chrom, aes(x=chrom$RTmins,y=chrom[,i]), colour = getCols[i], size = 0.25) 
    }
  }
  
  p1 <- p1 +  xlab('RT (mins)') + ylab('intensity')
  p1 <- p1 + theme(plot.title = element_text(hjust = 0.5, size = 18))
  p1 <- p1 + theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour="light grey", size=0.01),
                   panel.border = element_blank(),
                   panel.grid.major = element_line(colour="light grey", size=0.01),
                   panel.background = element_blank())
  p1 <- p1 + scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(chrom$RTmins), by = 1))
  p1 <- p1 + scale_y_continuous(expand = c(0, 0))
  p1
  
  ## save plot as png file
  ggsave(paste(fileName,".png"), plot = p1, device = NULL, path = NULL, scale = 1, width = 35, height = 22, dpi = 500,units = c("cm")) 
  
  ## if zoom is required
  if (zoom) {
    zoomPlot <- chrom[which(chrom$RTmins > zoomSection[1] && chrom$RTmins < zoomSection[2]),]
    
    zoomPlot <- ggplot(chromatogram, aes(x=RT, y=intensity)) 
    zoomPlot <- zoomPlot + geom_line(alpha=0.75, colour = "red", size = 0.2)
    
    zoomPlot <- zoomPlot + theme(plot.title = element_text(hjust = 0.5, size = 16), 
                                 axis.text=element_text(size=6), 
                                 axis.text.y=element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank(),
                                 axis.ticks.y=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 plot.margin = unit(c(1, 0, -0.5, -1), "lines"))  # margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    zoomPlot <- zoomPlot + xlab("") 
    zoomPlot <- zoomPlot + ylab("")
    zoomPlot <- zoomPlot + scale_y_continuous(expand = c(0, 0))
    #zoomPlot <- zoomPlot + theme(plot.margin=grid::unit(c(0,0,0,0), "null"))
    plot(zoomPlot) 
    
    ggsave(paste(fileName,"_ZOOM.png"), plot = zoomPlot, device = NULL, path = NULL, scale = 1, width = 12, height = 10, dpi = 350,units = c("cm"))
    
  }
}

