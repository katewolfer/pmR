pmRplotSpectrum <- function(chrom){

  ## round the values
  chrom$mz <- round(chrom$mz, 5)

  ## find 6 largest peaks
  chromCheck <- chrom[with(chrom,order(-intensity)),]
  chromCheck <- chromCheck[1:6,]
  sortChrom <- order(chrom$intensity, decreasing = TRUE)[1:6]

  chromPlot <- ggplot(chrom, aes(x = mz, y=intensity))
  #chromPlot <- chromPlot + geom_point(alpha=0.75, colour = "dark red", size = 0.2)
  chromPlot <- chromPlot + geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity),
                                        alpha=0.75, colour = "dark red")

  chromPlot <- chromPlot + geom_text(data=chrom[which(is.element(chrom$mz, chromCheck$mz)),],
                                     aes(label=mz, y=intensity), size=3, angle = 90,
                                     hjust = -0.25)


  chromPlot <- chromPlot + xlab("mass-to-charge (m/z)")
  chromPlot <- chromPlot + ylab("absolute intensity")
  chromPlot <- chromPlot + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.minor = element_line(colour="light grey", size=0.01),
                                 panel.border = element_blank(),
                                 panel.grid.major = element_line(colour="light grey", size=0.01),
                                 panel.background = element_blank())
  #chromPlot <- chromPlot + scale_x_continuous(breaks = round(seq(min(chromatogram$mz), max(chromatogram$mz), by = 1),1), expand = c(0, 0))
  chromPlot <- chromPlot + scale_y_continuous(expand = expansion(mult = .075))
  #chromPlot <- chromPlot + expand_limits(y = 0)
  #chromPlot <- chromPlot + scale_x_continuous(expand = c(0,0))
  chromPlot <- chromPlot + geom_hline(yintercept = 0)
  plot(chromPlot)
  #print(chrom)

  return(chromPlot)
}
