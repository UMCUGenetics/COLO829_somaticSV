#!/usr/bin/env Rscript


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggsci))

###PARSE ARGUMENTS####
option_list = list(
  make_option(c("-t", "--truthset"), type="character", default=NULL, 
              help="VCF file to use as truthset"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="VCF files to compare (input as comma-separated)"),
  make_option(c("-o", "--outPrefix"), type="character", default="out", 
              help="Prefix to produce output files (.tsv and .pdf) [out]"),
  make_option(c("-g", "--genome"), type="character", default="hg19", 
              help="Genome reference [hg19]"),
  make_option(c("-c", "--chromosome"), type="character", default=".", 
              help="Restrict analysis to a single chromosome or \".\" for all chromosomes [.]"),
  make_option(c("--maxgap"), type="integer", default=1000, 
              help="Maximum window around breakpoint for overlap with truth set"),
  make_option(c("--sizemargin"), type="integer", default=1, 
              help="Event size margin for overlap with truth set"),
  make_option(c("--margintosize"), type="integer", default=1, 
              help="Modify margin depending on size of the event for overlap with truth set")
); 

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

if (is.null(opt$truthset)){
  optparse::print_help(opt_parser)
  stop("Truthset VCF file must be specified with --truthset", call.=FALSE)
}
if (is.null(opt$file)){
  optparse::print_help(opt_parser)
  stop("At least one VCF file must be specified with --file", call.=FALSE)
}




fileList <- as.list(strsplit(opt$file, ",")[[1]])
outTsv <- paste(opt$outPrefix, ".tsv", sep = "")
outPdf <- paste(opt$outPrefix, ".pdf", sep = "")

truthsetFile <- opt$truthset



vcf2bpgr <- function(filename, genome, chromosome) {
  vcf <- VariantAnnotation::readVcf(filename, genome)
  bpgr <- StructuralVariantAnnotation::breakpointRanges(vcf, nominalPosition=T)
  if (chromosome != "."){
    bpgr <- bpgr[seqnames(bpgr) == chromosome |
                   seqnames(partner(bpgr, selfPartnerSingleBreakends=TRUE)) == chromosome]
  }
  names(bpgr) <- paste(filename, names(bpgr), sep = "_")
  bpgr$partner <- paste(filename, bpgr$partner, sep = "_")
  bpgr$caller <- filename
  return(bpgr)
}


# ###READ TRUTHSET###
# 
truthset_bpgr <- vcf2bpgr(filename=truthsetFile, genome=opt$genome, 
                          chromosome=opt$chromosome)

bpgrList <- sapply(fileList, vcf2bpgr, opt$genome, opt$chromosome)

svgr <- do.call("c", bpgrList)

svgr$truth_matches <- countBreakpointOverlaps(svgr, truthset_bpgr,
                                              maxgap = opt$maxgap,
                                              sizemargin=opt$sizemargin,
                                              restrictMarginToSizeMultiple = opt$margintosize,
                                              countOnlyBest = F,
                                              )

out.df <- as.data.frame(svgr) %>%
  dplyr::select(caller, truth_matches) %>%
  dplyr::group_by(caller) %>%
  dplyr::summarise(
    calls=dplyr::n(),
    tp=sum(truth_matches > 0),
  fp = calls - tp,
  fn = length(truthset_bpgr) - tp ) %>%
  dplyr::group_by(caller) %>%
  dplyr::mutate(
    Precision=tp / calls,
    Recall= tp / length(truthset_bpgr))
write.table(x = as.data.frame(out.df), file = outTsv, quote=F, sep="\t", row.names = F)

p <- ggplot(out.df) +
  aes(x=Recall, y=Precision, colour=caller) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=3) +
  theme_bw() +
  scale_color_npg(name="VCF file") +
  theme(text = element_text(face="bold", size = 12))
ggsave(filename = outPdf, plot = p, device = "pdf", width=10, height=8)
