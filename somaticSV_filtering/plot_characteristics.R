library(ggplot2)
library(ggrepel)

df <- read.table("truthset_characteristic_matrix.tsv",sep="\t")
colnames(df) <- c("tech", "type", "size", "value")
df$tech <- factor(df$tech, levels=c('ALL','ILLUMINA','NANOPORE','PACBIO','TENX'))
df$type <- factor(df$type, levels=c('DEL','INS','DUP','INV','BND'))
df$size <- factor(df$size, levels=c('<=100','100-1000','1000-10000','>10000'))


cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2")

p <- ggplot(data = df[df$tech == "ALL",], aes(x=size, y=value, fill = size, label=value)) +
  geom_bar(stat="identity", position = "dodge", col = "black") +
  geom_label(aes(y=value+1), position = position_dodge(width = .9),
             show.legend = F, fontface = "bold", col = "black") +
  facet_grid(.~type, space = "free_x", scales = "free_x", switch = "y") +
  
  scale_fill_manual(values = cbbPalette, 
                    limits =c('<=100','100-1000','1000-10000','>10000')) +
  xlab("") + ylab("# of SVs") + 
  theme_bw() +
  guides(fill=guide_legend(ncol=2)) + 
  scale_x_discrete(breaks=levels(df$size)) +
  scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5), position = "right") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 14, face="bold"),
        legend.key = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(colour=NA, fill=NA) 
  )

ggsave(filename = "truthset_characteristics.pdf", device="pdf", 
    plot = p, width = 8, height=5)
