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
packages <- c('tidyverse')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# test data for go_bubbl_plot.R
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
