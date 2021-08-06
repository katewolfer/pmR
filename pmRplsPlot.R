#df <- read.csv("AA summer PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("EPR summer PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("DTT summer PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("DCFH summer PLSR.csv", stringsAsFactors = FALSE)

#df <- read.csv("AA winter PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("EPR winter PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("DTT winter PLSR.csv", stringsAsFactors = FALSE)
#df <- read.csv("DCFH winter PLSR.csv", stringsAsFactors = FALSE)

getFiles <- list.files(pattern = "PLSR.csv")
getFiles <- getFiles[-1]

for (i in 1:length(getFiles)){
  df <- read.csv(getFiles[i],check.names = FALSE, stringsAsFactors = FALSE)
  df$date <- as.Date(df$date, format = "%d/%m/%Y")
  getSD1 <- sd(df$LV1, na.rm = TRUE)
  getSD2 <- getSD1*2
  getSD3 <- getSD1*3

  nb.cols <- 8
  mycolors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols))

  if(ncol(df) == 5){

    ## ONE latent variable
    plsPlot <- ggplot(df, aes(x = date, y = LV1, color = df[,4])) + geom_point(size = 5) + theme_bw()
    plsPlot <- plsPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plsPlot <- plsPlot + scale_color_gradientn(colours = mycolors)
    #plsPlot <- plsPlot + geom_ribbon(data = df, aes(x = runOrder, ymax=getSD1, ymin=0), fill="pink", alpha=.5)
    #plsPlot <- plsPlot + scale_x_discrete(labels = as.character(df$id))
    plsPlot <- plsPlot + scale_x_date(date_breaks = "2 day")
    plsPlot <- plsPlot + geom_hline(yintercept = 0)# + geom_text(aes(0,0,label = "2*SD", vjust = -1))
    plsPlot <- plsPlot + geom_hline(yintercept = getSD2, color = "orange", linetype = "dashed", size = 1)
    plsPlot <- plsPlot + geom_hline(yintercept = -getSD2, color = "orange", linetype = "dashed", size = 1)
    plsPlot <- plsPlot + geom_hline(yintercept = getSD3, color = "red", size = 1)
    #plsPlot <- plsPlot + geom_hline(yintercept = -getSD3, color = "red", size = 1)
    plsPlot <- plsPlot + labs(x = "date", y = "latent variable 1")
    plsPlot <- plsPlot + theme(text = element_text(size = 16))
    plsPlot <- plsPlot + theme(legend.position = "bottom",
                               legend.key.height = unit(0.3, "cm"),
                               legend.key.width = unit(2, "cm"),
                               legend.title = element_text(size = 16),
                               legend.text=element_text(size=14)) + theme(legend.title=element_blank())

  } else {

    ## TWO latent variables
    getHT2 <- as.data.frame(pcaMethods:::simpleEllipse(df$LV1, df$LV2, alfa = 0.95, len = 500))
    plsPlot <- ggplot(df, aes(x = LV1, y = LV2, color = df[,4])) + geom_point(size = 5) + theme_bw()
    plsPlot <- plsPlot + scale_color_gradientn(colours = mycolors)
    plsPlot <- plsPlot + labs(x = "latent variable 1", y = "latent variable 2")
    plsPlot <- plsPlot + theme(text = element_text(size = 16))
    plsPlot <- plsPlot + geom_point(data = getHT2, aes(x = V1, y = V2), colour = '#888888', size = 0.2)
    plsPlot <- plsPlot + geom_hline(yintercept = 0)
    plsPlot <- plsPlot + geom_vline(xintercept = 0)
    plsPlot <- plsPlot + theme(legend.position = "bottom",
                               legend.key.height = unit(0.3, "cm"),
                               legend.key.width = unit(2, "cm"),
                               legend.title = element_text(size = 16),
                               legend.text=element_text(size=14)) + theme(legend.title=element_blank())

  }

  pasteName <- gsub('PLSR.csv','PLSR plot.png', getFiles[i])
  ggsave(pasteName, plot = plsPlot, device = NULL, path = NULL,
         scale = 1, width = 16, height = 16, dpi = 200,units = c("cm"))

}
