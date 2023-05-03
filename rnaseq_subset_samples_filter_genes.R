#!/usr/bin/env Rscript

library('optparse')

# filter_functions <- c('mean', 'median', 'sd', 'sum')
# filter_tests <- c('lt', 'lte', 'gt', 'gte', 'eq')

option_list <- list(
  make_option("--output_file", type="character", default='all-subset.tsv',
              help="Output file name [default %default]" ),
  make_option("--samples_file", type="character", default=NULL,
              help="File of samples to subset to [default %default]" ),
  make_option("--genes_file", type="character", default=NULL,
              help="File of genes to subset to [default %default]" ),
  make_option("--no_counts", action="store_true", default=FALSE,
              help="Whether to exclude raw counts [default %default]" ),
  make_option("--no_normalised_counts", action="store_true", default=FALSE,
              help="Whether to exclude normalised counts [default %default]" ),
  # make_option("--filter_function", type="character", default=NULL,
  #             help=paste0("Measure to use to filter genes: ", 
  #                         paste0(filter_functions, collapse = ", "), " [default %default]" ) ),
  # make_option("--filter_test", type="character", default=NULL,
  #             help=paste0("Logical test to use to filter genes: ", 
  #                         paste0(filter_tests, collapse = ", "),  " [default %default]" ) ),
  # make_option("--filter_threshold", type="character", default=NULL,
  #             help="Threshold to use to filter genes [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to subset an all file to specific columns (Samples) and specific rows (Genes)',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'rnaseq_subset_samples_filter_genes.R',
    usage = "Usage: %prog [options] all_file",
    description = desc ),
  positional_arguments = 1
)

# cmd_line_args <-
#   list(options = list('output_file' = 'test_data/test_rnaseq_data.out.tsv',
#                       'samples_file' = 'test_data/test_control_samples.tsv',
#                       'genes_file' = 'test_data/test_genes.txt',
#                       'no_counts' = FALSE,
#                       'no_normalised_counts' = FALSE,
#                       # 'filter_function' = 'mean',
#                       # 'filter_test' = 'gt',
#                       # 'filter_threshold' = 1,
#                       debug = TRUE),
#        args = c('test_data/test_rnaseq_data.tsv'))

# check options
# if (!is.null(cmd_line_args$options[['filter_function']])) {
#   if (!(cmd_line_args$options[['filter_function']] %in% filter_functions)) {
#     stop('--filter_function must be one of ', paste0(filter_functions, collapse = ", ") )
#   }
# }
# if (!is.null(cmd_line_args$options[['filter_test']])) {
#   if (!(cmd_line_args$options[['filter_test']] %in% filter_tests)) {
#     stop('--filter_test must be one of ', paste0(filter_tests, collapse = ", ") )
#   }
# }

# load packages
packages <- c('tidyverse', 'rnaseqtools')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# unpack options
return_counts <- !cmd_line_args$options[['no_counts']]
return_normalised_counts <- !cmd_line_args$options[['no_normalised_counts']]

# load RNAseq data
data <- load_rnaseq_data(data_file = cmd_line_args$args[1])

# load samples file
samples <- load_rnaseq_samples(cmd_line_args$options[['samples_file']])

# subset to samples
subset_data <- 
  subset_to_samples(data, samples,
                    counts = return_counts,
                    normalised_counts = return_normalised_counts)

# # subset to genes
# norm_counts <- get_counts(subset_data, samples, normalised = TRUE)
# rowwise_means <- apply(norm_counts, 1, mean)
# subset_data <- filter(subset_data, rowwise_means > 1)

# subset to genes
if (!is.null(cmd_line_args$options[['genes_file']])) {
  genes <- read_tsv(cmd_line_args$options[['genes_file']])
  if (grepl("^ENS", colnames(genes)[1])) {
    genes <- read_delim(cmd_line_args$options[['genes_file']],
                        header = FALSE)
    colnames(genes)[1] <- "GeneID"
  }
  # join gene id column to subset_data
  filtered_data <- semi_join(subset_data, genes)
} else {
  filtered_data <- subset_data
}

write_tsv(filtered_data, file = cmd_line_args$options[['output_file']])
