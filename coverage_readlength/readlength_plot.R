library(ggplot2)
library(reshape2)
library(scales)
library(ggrepel)

#setwd("/home/cog/jvalleinclan/hpcMount/projects/TP0001_General/COLO/analysis/jvalleinclan/shared_COLO829/readlength/new")
#Make histograms!
mat <- matrix(ncol=4)

for (file in list.files(path="./", pattern="*readlength.hist$", recursive=TRUE) ) {
    print(file)
    data <- read.table(file, header=F, stringsAsFactors = F, sep=",")
    mat <- rbind(mat, data)
}
colnames(mat) <- c("bin","value","tech","type")
mat[mat == 0] <- NA
mat <- na.omit(mat)

values <- apply(mat, 1, function(x){
    tmp <- rep(x["bin"], times=round(as.numeric(x["value"])/1000))
})
values <- unlist(values, use.names = FALSE)
values <- as.numeric(values)


techs <- apply(mat, 1, function(x){
    tmp <- rep(x["tech"], times=round(as.numeric(x["value"])/1000))
})
techs <- unlist(techs, use.names = FALSE)


types <- apply(mat, 1, function(x){
    tmp <- rep(x["type"], times=round(as.numeric(x["value"])/1000))
})
types <- unlist(types, use.names = FALSE)

dist.df <- data.frame(tech=techs, value=values,type=types)
dist.df$type <- factor(dist.df$type, levels=c('normal', 'tumor'))
dist.df[dist.df$value == 0, "value"] <- 1


palette <- c('#d7191c','#fdae61','#ffffbf','#80cdc1','#2c7bb6')

mean_matrix <- aggregate(dist.df$value, by=list(dist.df$tech, dist.df$type), mean)
colnames(mean_matrix) <- c("tech", "type", "mean")

plot <- ggplot(dist.df, aes(x=value, fill=tech)) + 
    geom_density(alpha=.6) + 
    geom_label_repel(data=mean_matrix[mean_matrix$tech != "illumina",], aes(label=paste(round(mean/1000,0), "Kbp", sep = " "), x = mean, y = 4), 
                     col="black", show.legend = F, fontface = "bold", size = 6) +
  geom_label_repel(data=mean_matrix[mean_matrix$tech == "illumina",], aes(label=paste(round(mean,0), "bp", sep = " "), x = mean, y = 4), 
                   col="black", show.legend = F, fontface = "bold", size = 6) +
    facet_grid(type ~ .)  + 
    scale_fill_manual(values=palette) + 
    scale_color_manual(values=palette, guide=F) + 
    scale_x_log10(labels=comma, breaks = c(150,500, 1000,10000,25000,100000,500000), limits = c(100,1000000)) +
    ylim(c(0,6.5)) +
    xlab("Molecular length (bp)") + ggtitle("Single molecule lengths") +
  theme_bw() +
    theme(text = element_text(face="bold"),
          strip.background = element_rect(fill="white"),
          strip.text = element_text(size= 14, face="bold"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=14, face = "bold"),
          legend.text = element_text(size=14, face = "bold"),
          plot.title = element_text(size=14, face = "bold"),
          axis.title = element_text(size=14, face = "bold"))

plot
pdf("readlength_final.pdf", width = 10, height = 8)
plot
dev.off()

write.csv2(file="readlength_means.tsv", x = mean_matrix)
