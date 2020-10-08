library(ggplot2)
library(ggsci)

df <- read.table("NYGC_comparison_result.tsv", header = F, sep = "\t")
colnames(df) <- c("sequencer","result","count")
p <- ggplot(df, aes(x=sequencer, y=count, fill = result)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_npg() +
  scale_y_continuous(breaks=seq(0,80,10)) +
  theme_bw() + 
  xlab("") + ylab("") + labs(fill="") + 
  theme(
    text = element_text(size=14, face = "bold")
  )

ggsave(filename = "nygc_comparisons.pdf", plot = p, device = "pdf")
