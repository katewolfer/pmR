##################################################
## Tidy spectrum plotting function              ##
## Kate Wolfer, Institute of Hepatology, London ##
## 06 December 2018                             ##
##################################################

## INPUT
## chrom: csv file, with two columns: mz, showing mass, and
## intensity, showing intensity values for each mass
## cutStart: time in mins before which data should be discarded
## cutEnd: time in mins after which data should be discarded
## zoom: logical, 1 if subsection is required
## zoomSection: vector (e.g c(1,5)) indicating start and end of 
## RTmins to output zoom section of chromatogram

## OUTPUT

plotSpectrum <- function(chrom, cutStart, cutEnd, zoom, zoomSection, fileName) {
  require(ggplot2)
  require(RColorBrewer)
  
  chrom <- chrom[which(chrom$RTmins < cutStart),]
  chrom <- chrom[which(chrom$RTmins > cutEnd),]
  
  chromPlot <- ggplot(chrom, aes(x=mz, y=intensity)) 
  chromPlot <- chromPlot + geom_line(alpha=0.75, colour = "dark red", size = 0.75)
  chromPlot <- chromPlot + xlab("mass-to-charge ratio (m/z)")
  chromPlot <- chromPlot + ylab("intensity")
  chromPlot <- chromPlot + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.minor = element_line(colour="light grey", size=0.01),
                                 panel.border = element_blank(),
                                 panel.grid.major = element_line(colour="light grey", size=0.01),
                                 panel.background = element_blank())
  #chromPlot <- chromPlot + scale_x_continuous(breaks = round(seq(min(chromatogram$mz), max(chromatogram$mz), by = 1),1), expand = c(0, 0))
  chromPlot <- chromPlot + scale_y_continuous(expand = c(0, 0))
  chromPlot <- chromPlot + scale_x_continuous(expand = c(0, 0))
  plot(chromPlot) 
  
  ggsave(paste(fileName,".png"), plot = chromPlot, device = NULL, path = NULL, scale = 1, width = 17.5, height = 11, dpi = 500,units = c("cm"))
  
}