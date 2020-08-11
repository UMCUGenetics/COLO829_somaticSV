library(UpSetR)
filename <- "COLO829.somatic.overlap.somaticData.validations.upset.tsv"

df <- read.table(filename, sep="\t", header=T)
colnames(df) <- c("ILL", "PB", "10X","ONT", "PCR", "CAPTURE", "BNG", "MULT")

myfunction <- function(row, set1, set2, set3){
    newData <- (any(c(row[set1], row[set2], row[set3]) == 1))
}
pdf(file=paste(filename,".pdf", sep = ""), width = 16, height = 10)
upset(df, sets = rev(c("ILL", "ONT", "PB", "10X")), keep.order = T,
      point.size = 4, line.size=1.5, mainbar.y.label = "# Breakpoints",
      sets.x.label = "# of breakpoints", text.scale = 3.5,
     order.by = "freq", set_size.show = T, set_size.scale_max=100)

upset(df, sets = rev(c("PCR", "CAPTURE", "BNG", "MULT")), keep.order = T,
      # queries = list(list(query = myfunction, params = list("PCR", "CAPTURE", "BN"), color = "darkblue", active = T, query.name = "Validated")),
      # query.legend = "top",
      point.size = 4, line.size=1.5, mainbar.y.label = "# Breakpoints",
      sets.x.label = "# of breakpoints", text.scale = 3.5,
      order.by = "freq", set_size.show = T, set_size.scale_max=100)
dev.off()
