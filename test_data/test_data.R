#!/usr/bin/env Rscript

library('optparse')
option_list <- list()

desc <- paste(
  '\nScript to create test data for the scripts in the repo',
  'Files will be created in the current working directory',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'test_data.R',
    usage = "Usage: %prog [options] ",
    description = desc ),
  positional_arguments = 0
)

# load packages
packages <- c('tidyverse', 'DESeq2')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# test data for go_bubble_plot.R
set.seed(7531)
num_terms <- 150
go_bubble_plot_test_data <- tibble(
  GO.ID = sprintf('GO:%07d', seq_len(num_terms)),
  Term = sprintf('Term%d', seq_len(num_terms)),
  Significant = sample(1:50, num_terms, replace = TRUE),
  pval = 10^-(rnorm(num_terms, mean = 5, sd = 1)),
  cat = sample(c('BP', 'CC', 'MF'), num_terms, replace = TRUE)
)
for (domain in c('BP', 'CC', 'MF')) {
  data <- go_bubble_plot_test_data %>% 
    filter(., cat == domain) %>% 
    select(., -cat)
  write_tsv(data, path = paste0(domain, '.sig.tsv'))
}

# test data for graph_counts_by_group_facet.R
set.seed(583)
tmp <- makeExampleDESeqDataSet(n = 10, m = 12, betaSD = 2)
tmp2 <- makeExampleDESeqDataSet(n = 10, m = 12, betaSD = 2, interceptMean = 10)
# adjust sample names
rownames(colData(tmp2)) <- paste0('sample', 13:24)
# create samples file
df <- rbind(colData(tmp), colData(tmp2))
set.seed(714)
samples_df <- data.frame(
  sample = rownames(df),
  condition = rep( rep(c('wt', 'mut'), each = 6), 2 ),
  treatment = rep(c('control', 'treated'), each = 12),
  sex = sample(c("F", "M"), 24, replace = TRUE)
)
write_tsv(samples_df, path = 'test_samples.tsv')

# pivot samples file to use as metadata
samples_df_long <- samples_df %>% 
  pivot_longer(., cols = condition:sex, names_to = "category",
               values_to = "value")
write_tsv(samples_df_long, path = 'test_samples_long.tsv')

tmp <- DESeq(tmp)
res <- results(tmp)
tmp2 <- DESeq(tmp2)
# create gene info
num_rows <- 10
set.seed(208)
starts <- sample(1:10000, num_rows)
test_all_data <- tibble(
  'Gene ID' = sprintf('ZFG%03d', seq_len(num_rows)),
  'chr' = sample(1:25, num_rows, replace = TRUE),
  'start' = starts,
  'end' = as.integer(starts + 100),
  'strand' = sample(c('1', '-1'), num_rows, replace = TRUE),
  'p value' = res$pvalue,
  'Adjusted p value' = res$padj,
  'Gene name' = paste0('gene-', seq_len(num_rows))
)
test_all_data$chr <- factor(test_all_data$chr, levels = unique(test_all_data$chr))
test_all_data$strand <- factor(test_all_data$strand)

# create counts file
counts_1 <- counts(tmp)
counts_2 <- counts(tmp2)
colnames(counts_1) <- sub("$", " count", colnames(counts_1))
colnames(counts_2) <- sub("$", " count", colnames(counts_2))
norm_counts_1 <- counts(tmp, normalized = TRUE)
norm_counts_2 <- counts(tmp2, normalized = TRUE)
colnames(norm_counts_1) <- sub("$", " normalised count", colnames(norm_counts_1))
colnames(norm_counts_2) <- sub("$", " normalised count", colnames(norm_counts_2))
test_rnaseq_data <- cbind(
  test_all_data,
  counts_1,
  counts_2,
  norm_counts_1,
  norm_counts_2
)
write_tsv(test_rnaseq_data, path = 'test_rnaseq_data.tsv')

# test data for GO barchart
num_terms <- 30
set.seed(682)
go_bubble_plot_test_data <- tibble(
  GO.ID = sprintf('GO:%07d', seq_len(num_terms)),
  Term = sprintf('Term%d', seq_len(num_terms)),
  FE = sample(seq(3,20), num_terms, replace = TRUE),
  log10p = rnorm(num_terms, mean = 5, sd = 1),
  pval = 10^-log10p,
  Category = sample(c('BP', 'CC', 'MF'), num_terms, replace = TRUE),
  Set = sample(c('Expt1', 'Expt2', 'Expt3'), num_terms, replace = TRUE),
  up_down = sample(c('Up', 'Down'), num_terms, replace = TRUE),
)
write_tsv(go_bubble_plot_test_data, path = 'test_data_go.tsv')

# test data for gsea_to_genes
set.seed(241)
gsea_report <- tibble(
  NAME = sprintf('Term%03d', 1:5),
  `GS<br> follow link to MSigDB` = sprintf('Term%03d', 1:5),
  `GS DETAILS` = "Details...",
  SIZE = sample(100:200, 5),
  ES = runif(5, min = 0.5, max = 1.3),
  NES = ES * 1.8,
  `NOM p-val` = 0,
  `FDR q-val` = c(0, 0.0002, 0.001, 0.049, 0.06),
  `FWER p-val` = 0,
  `RANK AT MAX` = sample(500:2000, 5),
  `LEADING EDGE` = "notes"
)
write_tsv(gsea_report, path = 'test_gsea_report.xls')

for (term in gsea_report$NAME) {
  gene_info <- tibble(
    NAME = sprintf('row_%d', 0:19),
    PROBE = sprintf('ZFG%03d', 0:19),
    `GENE SYMBOL` = sprintf('ZFG%03d', 0:19),
    GENE_TITLE = sprintf('gene%d | description', 0:19),
    `RANK IN GENE LIST` = sample(3:100, 20, replace = TRUE),
    `RANK METRIC SCORE` = rnorm(20, mean = 20, sd = 2),
    `RUNNING ES` = runif(20),
    `CORE ENRICHMENT` = rep(c('Yes', 'No'), each = 10)
  )
  write_tsv(gene_info, path = paste0(term, '.xls'))
}