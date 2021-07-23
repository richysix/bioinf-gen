suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(DESeq2))

options(readr.show_progress=FALSE)

library('optparse')

option_list <- list(
  make_option("--annotation_file", type="character", default='../annotation/annotation.txt',
              help="File containing gene annotation [default %default]" ),
  make_option("--counts_file", type="character", default='counts.txt',
              help="File containing raw counts [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to run DESeq2 from a counts file',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'deseq2-from-counts.R',
    usage = "Usage: %prog [options]",
    description = desc ),
  positional_arguments = 0
)

# Get samples
samples <- read_tsv("samples.txt")
names(samples)[1] <- "sample"
groups <- setdiff(colnames(samples), c("sample", "condition"))
if (length(groups) > 0) {
  for (col in c('condition', groups)) {
    samples[[col]] <- factor(samples[[col]], levels = unique(samples[[col]]))
  }
}

# Get comparison
comparison <- fromJSON("comparison.json")
samples_deseq2 <- samples
samples_deseq2$condition <- fct_collapse(
  samples_deseq2$condition,
  exp=comparison$experimental_conditions,
  con=comparison$control_conditions
)
for (col in c('condition', groups)) {
  samples_deseq2[[col]] <- fct_relabel(
    samples_deseq2[[col]],
    ~ str_replace_all(.x, "-", "_")
  )
}

# Get annotation
annotation <- read_tsv(
  cmd_line_args$options[['annotation_file']],
  col_names=c(
    "Gene", "Chr", "Start", "End", "Strand",
    "Biotype", "Name", "Description"
  ),
  col_types=readr::cols(
    Gene=readr::col_character(),
    Chr=readr::col_factor(),
    Start=readr::col_integer(),
    End=readr::col_integer(),
    Strand=readr::col_integer(),
    Biotype=readr::col_factor(),
    Name=readr::col_character(),
    Description=readr::col_character()
  )
)

# Get counts
counts <- read_tsv(cmd_line_args$options[['counts_file']])
names(counts)[1] <- "gene"
# sort count columns by samples df
counts <- counts %>% select(., gene, all_of(samples_deseq2$sample))

# DESeq2
if (length(groups) > 0) {
  design <- as.formula(paste("~", paste(groups, collapse = " + "), "+ condition"))
} else {
  design <- formula(~ condition)
}

dds <- DESeqDataSetFromMatrix(
  column_to_rownames(counts, var="gene"),
  column_to_rownames(samples_deseq2, var="sample"),
  design=design
)
rm(counts)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "exp", "con"), alpha=0.05, tidy=TRUE)

# Output
write_tsv(
  enframe(dds$sizeFactor),
  "size-factors.tsv",
  col_names=FALSE
)
res <- inner_join(
  select(res, Gene=row, pval=pvalue, adjp=padj, log2fc=log2FoldChange),
  annotation,
  by="Gene"
)
res <- inner_join(
  res,
  rownames_to_column(
    rename_all(
      as.data.frame(counts(dds)),
      paste0,
      " count"
    ),
    var="Gene"
  ),
  by="Gene"
)
res <- inner_join(
  res,
  rownames_to_column(
    rename_all(
      as.data.frame(counts(dds, normalized=TRUE)),
      paste0,
      " normalised count"
    ),
    var="Gene"
  ),
  by="Gene"
)
res <- arrange(res, adjp, pval, Gene)
write_tsv(res, "all.tsv")
write_tsv(filter(res, adjp < 0.05), "sig.tsv")

# QC
#rld <- rlog(dds, blind=TRUE)
rld <- vst(dds, blind=TRUE)
colData(rld)$condition <- samples$condition
pdf("qc.pdf")
if (length(groups) > 0) {
  plotPCA(rld, intgroup=groups)
} else {
  plotPCA(rld)
}
plotMA(dds)
plotDispEsts(dds)
dev.off()
