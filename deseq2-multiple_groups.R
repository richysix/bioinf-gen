suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(DESeq2))

options(readr.show_progress=FALSE)

# Get samples
samples <- read_tsv("samples.txt")
groups <- setdiff(colnames(samples), c("sample", "condition"))
for (col in c('condition', groups)) {
  samples[[col]] <- factor(samples[[col]], levels = unique(samples[[col]]))
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
  "../annotation/annotation.txt",
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
counts <-tibble(
  gene=character(),
  sample=factor(samples$sample, levels=unique(samples$sample))[0],
  stranded2=integer(),
  .rows=0
)
for (sample in samples$sample) {
  tmp_counts <- read_tsv(
    str_c("../star2/", sample, "/ReadsPerGene.out.tab"),
    skip=4,
    col_names=c("gene", "unstranded", "stranded1", "stranded2"),
    col_types=cols(
      gene=col_character(),
      unstranded=col_skip(),
      stranded1=col_skip(),
      stranded2=col_integer()
    )
  )
  stop_for_problems(tmp_counts)
  tmp_counts <- add_column(
    tmp_counts,
    sample=factor(sample, levels=levels(counts$sample)),
    .before="stranded2"
  )
  counts <- bind_rows(counts, tmp_counts)
}
rm(tmp_counts)

# Aggregate counts
counts <- pivot_wider(counts, names_from=sample, values_from=stranded2)

# DESeq2
design <- as.formula(paste("~", paste(groups, collapse = " + "), "+ condition"))
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
plotPCA(rld, intgroup=groups)
plotMA(dds)
plotDispEsts(dds)
dev.off()
