library(ggplot2)
library(ggsci)
#Doing this through COVID-19, in shady coding conditions
sessionInfo()
df = read.table("final_somatic_list.tsv", header = T, sep ="\t", stringsAsFactors=F)
df$samplot <- factor(df$samplot, levels = c("False_positive", "Germline",
                                            "Duplicated", "Somatic"))
p <- ggplot(df, aes(x=samplot, fill = samplot)) +
  geom_bar(stat="count", alpha = .85, col = "black") +
  geom_text(stat='count', aes(label=..count..), 
            vjust=1.5,  fontface = "bold", col = "black", size = 5)+
  theme_bw() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10), expand = c(0,0)) + 
  scale_fill_npg() + xlab("") + ylab("") +
  theme(
    text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

pdf("samplot_results.pdf")
p
dev.off()
