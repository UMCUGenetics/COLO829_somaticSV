library(UpSetR)

filename <- "somatic.merged.clean.validations.upset.tsv"
df <- read.table(filename, sep="\t", header=T)
colnames(df) <- c("ILL", "PB", "10X","ONT", "PCR", "CAPTURE", "BNG")

df$MULT_NOVAL <- apply(df, 1, function(x){
    
    if (x["ILL"] + x["PB"] + x["ONT"] + x["10X"] > 1 & all(c(x["PCR"], x["CAPTURE"], x["BNG"]) == 0)) {
        res = 1
    } else {
        res = 0
    }
    res
})
df$SING_NOVAL <- apply(df, 1, function(x){
    if (x["ILL"] + x["PB"] + x["ONT"] + x["10X"] == 1 & all(c(x["PCR"], x["CAPTURE"], x["BNG"]) == 0)) {
        res = 1
    } else {
        res = 0
    }
    res
})

myfunction <- function(row, set1, set2, set3){
    newData <- (any(c(row[set1], row[set2], row[set3]) == 1))
}
pdf(file="somatic.merged.clean.upset.uniquified.pdf", width = 16, height = 10)
upset(df, sets = rev(c("ILL", "ONT", "PB", "10X")), keep.order = T,
      point.size = 4, line.size=1.5, mainbar.y.label = "# Breakpoints",
      sets.x.label = "# of breakpoints", text.scale = c(2, 2, 2, 2, 2, 2),
     order.by = "freq",  set_size.show = T, set_size.scale_max=7400)
upset(df, sets = rev(c("PCR", "CAPTURE", "BNG", "MULT_NOVAL")), keep.order = T,
      point.size = 4, line.size=1.5, mainbar.y.label = "# Breakpoints",
      sets.x.label = "# of breakpoints", text.scale = c(2, 2, 2, 2, 2, 2),
      order.by = "freq", set_size.show = T, set_size.scale_max=150)
dev.off()
