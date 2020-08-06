library(gdata)
library(reshape2)

nameFromPath <- function(path)
{
    name.list <- unlist(strsplit(path, '/'))
    name.tmp <- name.list[length(name.list)]
    name.tmp <- unlist(strsplit(name.tmp, '\\.'))[1]
    name <- unlist(strsplit(name.tmp, '_'))[2]
    name
}

files <- list.files(path=dir, pattern = "data_*", full.names = T) 
filenames <- lapply(files, nameFromPath)
names(files) <- filenames
df.list <- lapply(filenames, function(x) {
    tmp <- read.table(files[x], stringsAsFactors = F, sep = '\t', header = T)
    tmp$tech <- x
    tmp.melt <- melt(tmp, id.var = c("chrom", "start", "end", "sv_type", "tech"))
    tmp.melt <- tmp.melt[tmp.melt$value != 0 & tmp.melt$variable != "Total.Scores",]
    tmp.melt$value <- NULL
    tmp.melt
})
names(df.list) <- filenames





df <- Reduce(function(x, y) merge(x, y, by = c("chrom", "start", "end", "sv_type", "tech", "variable"), all = TRUE), df.list)
df <- data.frame(lapply(df, function(x){gsub("\\.\\.\\.\\.", "", x)}))
df <- reshape(df, idvar = c("chrom", "start", "end", "sv_type"), timevar = "tech", direction = "wide")
colnames(df) <- gsub("variable\\.", "", colnames(df))
write.table(df, file="merged_samplot_results.tsv", sep = '\t', quote =F, row.names = F)
