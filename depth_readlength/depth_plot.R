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
plot_matrix$tech <- as.character(plot_matrix$tech)
plot_matrix$tech[plot_matrix$tech == "10x"] <- "10X"
plot_matrix$tech[plot_matrix$tech == "illumina"] <- "ILL"
plot_matrix$tech[plot_matrix$tech == "bionano"] <- "BNG"
plot_matrix$tech[plot_matrix$tech == "ont"] <- "ONT"
plot_matrix$tech[plot_matrix$tech == "pacbio"] <- "PB"

plot_matrix$tech <- factor(plot_matrix$tech, levels=c("ILL", "ONT", "PB", "10X", "BNG"))

plot_matrix$type <- as.character(plot_matrix$type)
plot_matrix$type[plot_matrix$type == "tumor"] <- "COLO829"
plot_matrix$type[plot_matrix$type == "normal"] <- "COLO829BL"


mean_matrix <- aggregate(plot_matrix$coverage, by=list(plot_matrix$tech, plot_matrix$type), mean)
colnames(mean_matrix) <- c("tech", "type", "mean")

palette <- c('#2c7bb6', '#d7191c')

#Geomsplitviolifn from https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

plot <- ggplot(plot_matrix, aes(x=tech, fill=type, y=coverage)) + 
  geom_split_violin(draw_quantiles=c(.5), lwd = .8) + 
  #facet_grid(type ~ .) + 
  scale_y_continuous(limits = c(0,300), breaks = seq(0,300,50)) +
  scale_fill_manual(values=palette) + 
  xlab("Technology") + ylab("Depth") +
  theme_bw() +
  theme(text=element_text(size= 14, face="bold"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size= 14, face="bold"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14, face = "bold"),
        plot.title = element_text(size=14, face = "bold"),
        axis.title = element_text(size=14, face = "bold"),
        legend.title = element_blank())
  

pdf("depth_final.pdf", width = 10, height = 8)
plot
dev.off()
write.csv2(file="depth_means.tsv", x = mean_matrix)

