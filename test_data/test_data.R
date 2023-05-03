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
packages <- c('tidyverse', 'DESeq2', 'rprojroot')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

root_path <- find_root(is_rstudio_project)

# test data for GO kappa scores
num_sig_genes <- 100
sig <- tibble(
  Gene = paste0('gene', seq_len(num_sig_genes))
)
write_tsv(sig, file = file.path(root_path, 'test_data', 'sig.tsv'))
num_terms <- 10
sig_genes_per_term <- 50
set.seed(637)
go_sig_genes <- tibble(
  GO.ID = rep(sprintf('GO:%07d', seq_len(num_terms)), each = sig_genes_per_term),
  Term = rep(sprintf('Term%d', seq_len(num_terms)), each = sig_genes_per_term),
  Gene = sample(paste0('gene', seq_len(num_sig_genes)*2), num_terms*50, replace = TRUE)
) %>% 
  mutate(., `p value` = case_when(Gene %in% sig$Gene ~ 1,
                                  TRUE ~ 0)) %>% 
  unique()
write_tsv(go_sig_genes, file = file.path(root_path, 'test_data', 'BP.sig.genes.tsv'))

# make a file with just the control samples in
# load samples file
samples_df <- read_tsv(file.path(root_path, 'test_data', 'test_samples.tsv'))
samples_df %>%
  filter(treatment == "control") %>% 
  write_tsv(., file = file.path(root_path, 'test_data', 'test_samples_control.tsv'))

# pivot samples file to use as metadata
samples_df %>% 
  pivot_longer(., cols = condition:sex, names_to = "category",
               values_to = "value")
write_tsv(samples_df_long, file = file.path(root_path, 'test_data', 'test_samples_long.tsv'))

# load test_all_data
test_all_data <- read_tsv(file = file.path(root_path, 'test_data', 'test_rnaseq_data.tsv')) %>% 
  dplyr::select(GeneID:GO_MF)

# subset to 3 genes
test_all_data %>% filter(., GeneID %in% c('ENSTEST005', 'ENSTEST006', 'ENSTEST009')) %>% 
  select(., GeneID, Name = `Gene name`) %>% 
  write_tsv(., file = file.path(root_path, 'test_data', 'test_genes_to_label.txt'))

# make gene metadata file
gene_metadata <- select(test_all_data, GeneID, Class, GO_BP, GO_CC, GO_MF) %>% 
  pivot_longer(., -GeneID, names_to = 'category') %>% 
  filter(., !is.na(value)) %>% 
  arrange(., category, value)

write_tsv(gene_metadata, file = file.path(root_path, 'test_data', 'test_gene_metadata.tsv'))

set.seed(912)
sample_n(test_all_data, 3) %>% select(., GeneID) %>% 
  write_tsv(., file = file.path(root_path, 'test_data', 'test_genes.txt'))

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
write_tsv(go_bubble_plot_test_data, file = file.path(root_path, 'test_data', 'test_data_go.tsv'))

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
write_tsv(gsea_report, file = file.path(root_path, 'test_data', 'test_gsea_report.xls'))

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
  write_tsv(gene_info, file = file.path(root_path, 'test_data', paste0(term, '.xls')))
}
write.table(gene_info$PROBE[sample(1:20, 10)], quote = FALSE,
            file = file.path(root_path, 'test_data', 'gsea-genes.txt'), 
            col.names = FALSE, row.names = FALSE)

# make test data for upset-sig-genes.R
num_genes <- 100
genes <- sprintf('ENSTEST%03d', seq_len(num_genes))
names(genes) <- paste0('gene-', seq_len(num_genes))
set.seed(9374)
sets <- purrr::map(1:3, function(num){
  num_rows <- 20*num
  starts <- sample(1:10000, num_rows)
  gene_subset <- sample(genes, num_rows)
  tibble(
    'GeneID' = gene_subset,
    'chr' = sample(1:25, num_rows, replace = TRUE),
    'start' = starts,
    'end' = as.integer(starts + 100),
    'strand' = sample(c('1', '-1'), num_rows, replace = TRUE),
    'p value' = runif(num_rows, min = 1e-8, max = 0.001),
    'Adjusted p value' = `p value` * 20,
    'Gene name' = names(gene_subset)
  )
})
invisible(purrr::map(1:3, ~ write_tsv(sets[[.x]], file = paste0('set', .x, '.sig.tsv'))))
