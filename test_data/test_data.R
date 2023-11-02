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
