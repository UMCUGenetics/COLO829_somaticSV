library(ggplot2)
library(ggrepel)

plot_matrix <- matrix(ncol=3)[-1,]

for (file in list.files(path="./", pattern="*depth$", recursive=TRUE) ) {
  print(file)
  tech <- strsplit(file,"\\.")[[1]][1]
  type <- strsplit(file,"\\.")[[1]][2]
  print(tech)
  print(type)
  if (tech == 'bionano') {
    data <- read.table(paste(file),header=F)
    print(paste(tech, type, mean(data$V4), median(data$V4), sd(data$V4)))
    m <- cbind(type, cbind(tech, data$V4))
  } else {
    data <- read.table(file,header=T)
    print( paste( tech, type, mean(data$COV), median(data$COV), sd(data$COV)))
    m <- cbind(type, cbind(tech, data$COV))
  }
  
  plot_matrix <- rbind(plot_matrix, m)
}

colnames(plot_matrix)[3] <- "coverage"
plot_matrix <- as.data.frame(plot_matrix)
plot_matrix$coverage <- as.numeric(as.vector(plot_matrix$coverage))


mean_matrix <- aggregate(plot_matrix$coverage, by=list(plot_matrix$tech, plot_matrix$type), mean)
colnames(mean_matrix) <- c("tech", "type", "mean")

palette <- c('#d7191c','#fdae61','#ffffbf','#80cdc1','#2c7bb6')
plot <- ggplot(plot_matrix, aes(x=coverage,fill=tech)) + 
  geom_density(alpha=.6) +
  facet_grid(type~.) +
  geom_label_repel(data=mean_matrix, aes(label=paste(round(mean,0),"x", sep = " "), x = mean, y = 0.04), 
                   col="black", show.legend = F, fontface = "bold", size = 6) +
  scale_fill_manual(values=palette) +
  scale_x_continuous(limits=c(0,300)) +
  xlab("Coverage") + ggtitle("Coverage") +
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
          strip.text = element_text(size= 14, face="bold"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=14, face = "bold"),
          legend.text = element_text(size=14, face = "bold"),
          plot.title = element_text(size=14, face = "bold"),
          axis.title = element_text(size=14, face = "bold")) 

pdf("depth_final.pdf", width = 10, height = 8)
plot
dev.off()
write.csv2(file="depth_means.tsv", x = mean_matrix)

