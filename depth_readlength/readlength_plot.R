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
dist.df$tech <- as.character(dist.df$tech)
dist.df$tech[dist.df$tech == "tenx"] <- "10X"
dist.df$tech[dist.df$tech == "illumina"] <- "ILL"
dist.df$tech[dist.df$tech == "bionano"] <- "BNG"
dist.df$tech[dist.df$tech == "nanopore"] <- "ONT"
dist.df$tech[dist.df$tech == "pacbio"] <- "PB"

dist.df$tech <- factor(dist.df$tech, levels=c("ILL", "ONT", "PB", "10X", "BNG"))

dist.df$type <- as.character(dist.df$type)
dist.df$type[dist.df$type == "tumor"] <- "COLO829"
dist.df$type[dist.df$type == "normal"] <- "COLO829BL"
dist.df[dist.df$value == 0, "value"] <- 1


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
mean_matrix <- aggregate(dist.df$value, by=list(dist.df$tech, dist.df$type), mean)
colnames(mean_matrix) <- c("tech", "type", "mean")



plot <- ggplot(dist.df, aes(x=tech, fill=type, y=value)) + 
  geom_split_violin(draw_quantiles=c(.5), lwd=.8) + 
  #facet_grid(type ~ .) + 
  scale_y_log10(labels=comma, breaks = c(100,500, 1000, 5000, 10000, 50000,100000,500000), limits = c(100,1000000)) +
  scale_fill_manual(values=palette) + 
  xlab("Technology") + ylab("Molecular length (bp)") +
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

pdf("readlength_final.pdf", width = 10, height = 8)
plot
dev.off()

write.csv2(file="readlength_means.tsv", x = mean_matrix)
