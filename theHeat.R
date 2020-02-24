theHeat         <- function(inputMat,title=NA, xaxis, yaxis){
  p <- qplot(x=Var1, y=Var2, data=melt(cor(inputMat, use="pairwise.complete.obs")), fill=value, geom="tile") + scale_fill_gradient2(limits=c(-1, 1),low="darkblue",high="darkred") + ggtitle(title) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + theme(axis.title = element_blank())
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + theme(axis.text.x=element_text(hjust=0.95,vjust=0.2))
  return(p)  
}  