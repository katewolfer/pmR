## PCA version
library(pcaMethods)
library(ggforce)
df <- read.csv("Beijing PCA MASS.csv")
df <- read.csv("Beijing PCA MASS loadings.csv")

getHT2 <- as.data.frame(pcaMethods:::simpleEllipse(df$PC1, df$PC2, alfa = 0.95, len = 500))

nb.cols <- nrow(df)
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

mycolors <- df$season
mycolors[which(mycolors == "winter")] <- '#1338BE'
mycolors[which(mycolors == "summer")] <- '#F9A602'

seasonPlot <- ggplot(df, aes(x = PC1, y = PC2, color = Season))
seasonPlot <- seasonPlot + geom_point(size = 5)#stat = 'identity',colour = mycolors, size = 5)
seasonPlot <- seasonPlot + scale_colour_manual(values = c('#F9A602','#1338BE'))
seasonPlot <- seasonPlot + geom_point(data = getHT2, aes(x = V1, y = V2), colour = '#888888', size = 0.2)
seasonPlot <- seasonPlot + geom_hline(yintercept = 0) + labs(legend = "TITULO")
seasonPlot <- seasonPlot + geom_vline(xintercept = 0) + theme_bw()
seasonPlot <- seasonPlot + labs(x = "principal component 1", y = "principal component 2")
seasonPlot <- seasonPlot + theme(text = element_text(size = 20))
seasonPlot <- seasonPlot + theme(legend.position = "bottom") + theme(legend.title=element_blank())
seasonPlot

ggsave("beijing season PCA plot 2.png", plot = seasonPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))



nb.cols <- 15
mycolors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols))
plotlabel = bquote('AA (counts '*mu~g^-1*')')
getBreaks = seq(0,800,50)

plotlabel = bquote('EPR (counts '*mu~g^-1*')')
getBreaks = seq(min(df$EPR,na.rm = TRUE),max(df$EPR,na.rm = TRUE),1500)

plotlabel = bquote('DTT (nmol min'^-1 *mu~g^-1*')')
getBreaks = round(seq(min(df$DTT,na.rm = TRUE),max(df$DTT,na.rm = TRUE),5)+0.1,2)

plotlabel = bquote('DCFH (nmol ' ~H[2] ~O[2] *mu~g^-1*')')
getBreaks = round(seq(min(df$DCFH,na.rm = TRUE),max(df$DCFH,na.rm = TRUE),0.0025),4)


aaPlot <- ggplot(df, aes(x = PC1, y = PC2, color = DTT))
aaPlot <- aaPlot + geom_point(size = 5)
aaPlot <- aaPlot + scale_color_gradientn(colours = mycolors, breaks = getBreaks)
aaPlot <- aaPlot + geom_point(data = getHT2, aes(x = V1, y = V2), colour = '#888888', size = 0.2)
aaPlot <- aaPlot + geom_hline(yintercept = 0)
aaPlot <- aaPlot + geom_vline(xintercept = 0) + theme_bw()
aaPlot <- aaPlot + theme(text = element_text(size = 20))
aaPlot <- aaPlot + labs(colour=plotlabel)
aaPlot <- aaPlot + labs(x = "principal component 1", y = "principal component 2")
aaPlot <- aaPlot + theme(legend.position = "bottom",
                         legend.key.height = unit(0.3, "cm"),
                         legend.key.width = unit(4, "cm"),
                         legend.title = element_text(size = 16),
                         legend.text=element_text(size=14)) + theme(legend.title=element_blank())
#aaPlot <- aaPlot + guides(colour = guide_legend(title.position = "right"))
aaPlot

ggsave("beijing AA PCA plot 3.png", plot = aaPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))
ggsave("beijing EPR PCA plot 3.png", plot = aaPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))
ggsave("beijing DTT PCA plot 3.png", plot = aaPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))
ggsave("beijing DCFH PCA plot 3.png", plot = aaPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))

## all other PCAs
pcaPlot <- ggplot(df, aes(x = PC1, y = PC2, colour = EPR))
pcaPlot <- pcaPlot + geom_point(size = 5)#stat = 'identity',colour = mycolors, size = 5)
pcaPlot <- pcaPlot + scale_color_gradientn(colours = mycolors)
pcaPlot <- pcaPlot + geom_point(data = getHT2, aes(x = V1, y = V2), colour = '#888888', size = 0.2)
pcaPlot <- pcaPlot + geom_hline(yintercept = 0)
pcaPlot <- pcaPlot + geom_vline(xintercept = 0) + theme_bw()
pcaPlot <- pcaPlot + labs(x = "principal component 1", y = "principal component 2")
pcaPlot <- pcaPlot + theme(text = element_text(size = 20))
pcaPlot

ggsave("beijing EPR PCA plot.png", plot = pcaPlot, device = NULL, path = NULL,
       scale = 1, width = 26, height = 17, dpi = 600,units = c("cm"))

nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(11, "Dark2"))(nb.cols)
names(mycolors) <- levels(df$type)

mycolors <- c("#ffa500", # AMS
              "#b1f2ff", # biomass
              "#b11226", # cooking
              "#C0C0C0", # gas radical
              "#ffbeb1", # gases
              "#ffef00", # metal

              "#0dd9f4", # meteo
              "#6bafc3", # n-alkane
              "#ff721a", # PAH
              "#000000", # photolysis
              "#6091e3", # small ions

              "#008141", # SOA
              "#a2df35", # total EC
              "#ccb2e5", # total OC
              "#7f20b4", # vehicle
              "#f66f6f") # VOC

library(ggrepel)

## loadings plot
pcaPlot <- ggplot(df, aes(x = PC1, y = PC2, colour = type))
pcaPlot <- pcaPlot + geom_point(size = 5.5)
pcaPlot <- pcaPlot + scale_colour_manual(values = mycolors)
pcaPlot <- pcaPlot + geom_text_repel(aes(label=id),vjust=-2, color = "#888888", size = 3.5)
pcaPlot <- pcaPlot + geom_hline(yintercept = 0)
pcaPlot <- pcaPlot + geom_vline(xintercept = 0) + theme_bw()
pcaPlot <- pcaPlot + labs(x = "principal component 1", y = "principal component 2")
pcaPlot <- pcaPlot + theme(text = element_text(size = 20))
pcaPlot <- pcaPlot + theme(legend.title=element_blank())
pcaPlot <- pcaPlot + theme(legend.position="bottom")
pcaPlot

ggsave("beijing PCA loadings plot with labels.png", plot = pcaPlot, device = NULL, path = NULL,
       scale = 1, width = 29, height = 29, dpi = 500,units = c("cm"))
