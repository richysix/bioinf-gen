#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='histogram.pdf',
              help="Name of output file [default %default]" ),
  make_option("--data_column", type="character", default=NULL,
              help="Name of column to plot [default %default]" ),
  make_option("--binwidth", type="numeric", default=NULL,
              help="Name of column to plot [default %default]" )
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'histogram.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

# packages
packages <- c('tidyverse')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# unpack options
output_file <- cmd_line_args$options[['output_file']]
data_column <- cmd_line_args$options[['data_column']]
binwidth    <- cmd_line_args$options[['binwidth']]

# load data
plot_data <- read_tsv(cmd_line_args$args[1])

# try and find column
if (!is.null(data_column)) {
  if (!(hist_column %in% names(plot_data)) ) {
    stop('--data_column does not exist in the data\n')
  }
} else {
  # use first numeric column
  
}

# check column is numeric
if (class(plot_data[[hist_column]]) != "numeric") {
  stop('--data_column is not numeric\n')
}

# plot histogram
names(plot_data)[ names(plot_data) == data_column ] <- 'xvar'
histogram <- ggplot(data = plot_data, aes(x = xvar)) +
  geom_histogram(binwidth = binwidth)

# output plot
pdf(output_file)
print(histogram)
dev.off()

