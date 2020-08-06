library(corrplot)


data <- read.table("mergedCNA.tsv",sep="\t",header=T)
data <- data[,seq(4,8,1)]


pdf("correlation_plot.pdf", height =10, width = 10)
corrplot(cor(data),method='number',shade.col=NA,
         type = "upper", number.cex = 1.5, number.digits = 2,
         tl.col='black',tl.cex = 1.5,
         tl.srt=45,addCoef.col = 'black',insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(200),
         #cl.cex = 1.2, cl.lim = c(0,1), cl.align.text = "l",
         cl.pos="n"
         )
dev.off()
