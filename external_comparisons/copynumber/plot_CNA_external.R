library("ggplot2")
library("gtools")
library(reshape2)
#Largely inspired by https://github.com/UMCUGenetics/IAP/blob/develop/scripts/makeKaryotype.R
#from @Joep de Ligt and @Robert Ernst


myCols <- c("darkgray","steelblue3","red3","lightgray")
names(myCols) <- c("neutral1","gain","loss","neutral2")

df <- read.table("merged_external.tsv", header = T , sep = "\t", stringsAsFactors = F)

df$chr <- gsub(pattern = "hs", replacement = "", df$chr)

chromosomes <- mixedsort(as.character(unique((df$chr))))
dfm <- melt(df, id.vars=c("chr", "start", "end"))
dfm$value <- dfm$value+2
dfm$value[dfm$value > 6] <- 6
dfm$chr <- factor(dfm$chr, levels=chromosomes)
dfm$position <- (dfm$start + dfm$end)/2
dfm$cols <- apply(dfm, 1, function(x){
  if (as.integer(x["value"]) > 2 ) {
    col <- "gain"
  } else if (as.integer(x["value"]) < 2) {
    col <- "loss"
  } else {
    col <- 'neutral1'
  }
  #   if (x["chr"] == "X") {
  #     col <- 'neutral1'
  #   }
  #   else if (as.integer(x["chr"])%%2 == 0){
  #     col <- 'neutral1'
  #   } else {
  #     col <- 'neutral2'
  #   }
  # }
  col
})

p <- ggplot(dfm, aes(position, value, group=chr)) + 
  geom_point(aes(colour=cols), shape=20, size=2) + 
  scale_y_continuous(limits = c(0,6), breaks =seq(0,6,1)) +  
  facet_grid(variable~chr, scales="free_x", space = "free_x") + 
  scale_colour_manual(name="event", values=myCols) +
  xlab("") + ylab("Copy number") +
  theme(
    text=element_text(size=16, face="bold"),
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.2),
    
    legend.position="none",
    panel.background=element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background=element_blank(),
    strip.background = element_rect(fill=NA,colour="black", size=1),
    panel.spacing.x = unit(0, "lines"),
    panel.spacing.y = unit(0.2, "lines")
  )
p
ggsave(plot = p, filename = "cna_plot_external.pdf", device="pdf", width=14, height = 8)
