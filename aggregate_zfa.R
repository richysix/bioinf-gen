#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--group_by_cols", type="character", action="store", default='Set,Parent.Term',
              help="Name of the columns to group by for aggregation [default %default]"),
  make_option("--filter_by_cols", type="character", action="store", default='Parent.Term',
              help="Name of the columns to filter by [default %default]"),
  make_option("--log10p_col", type="character", action="store", default='log10p',
              help="Name of the column for log10(pvalue) [default %default]"),
  make_option("--fe_col", type="character", action="store", default='FE',
              help="Name of the column for Fold Enrichment [default %default]")
)

desc <- paste('', 'Script to aggregate ZFA terms by another column',
    '(e.g. Parent.Term)',
    sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'aggregate_zfa.R',
    usage = "Usage: %prog [options] input_file output_file",
    description = desc),
  positional_arguments = 2
)

# unpack options
group_by_cols <-
        unlist( strsplit(cmd_line_args$options[['group_by_cols']],
                          ',', fixed = TRUE) )
filter_by_cols <-
        unlist( strsplit(cmd_line_args$options[['filter_by_cols']],
                          ',', fixed = TRUE) )
log10p_col <- cmd_line_args$options[['log10p_col']]
fe_col <- cmd_line_args$options[['fe_col']]

packages <- c('tidyverse')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

data <- read.delim(cmd_line_args$args[1])

# set levels of group by columns
for (col in group_by_cols) {
    data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
}

# convert log10 and fe col names to expected col names
colnames(data)[ colnames(data) == log10p_col ] <- 'log10p'
colnames(data)[ colnames(data) == fe_col ] <- 'FE'

aggregated_data <- data %>%
    filter_at(filter_by_cols, all_vars(!is.na(.))) %>%
    group_by_at(group_by_cols) %>%
    summarise(max_log10p = max(log10p), max_fe = max(FE), numTerms = n()) %>%
    arrange(desc(numTerms), .by_group = TRUE)

write.table(aggregated_data, file = cmd_line_args$args[2], quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)
