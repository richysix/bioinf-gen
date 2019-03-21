#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  #make_option("--group_by_cols", type="character", action="store", default='Set,Parent.Term',
  #            help="Name of the column to group by for aggregation [default]")
  #make_option("--arrange_cols", type="character", action="store", default='Cluster',
  #            help="Name of the column to group by for aggregation [default]")
)

desc <- paste('', 'Script to aggregate ZFA terms by another column',
    '(e.g. Parent.Term)',
    paste('At the moment the script assumes the columns Set and Parent.Term',
           'exist to group by.'),
    'Eventually it will allow those columns to be specified.',
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
#group_by_cols <-
#    as.list(
#        unlist( strsplit(cmd_line_args$options[['group_by_cols']],
#                          ',', fixed = TRUE) )
#    )
#group_by_cols_q <- quos(group_by_cols)

packages <- c('tidyverse')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

data <- read.delim(cmd_line_args$args[1])
aggregated_data <- data %>%
    filter(!is.na(Parent.Term)) %>%
    group_by(Set, Parent.Term) %>%
    summarise(max_log10p = max(log10p), max_fe = max(FE), numTerms = n()) %>%
    arrange(desc(numTerms), .by_group = TRUE)

write.table(aggregated_data, file = cmd_line_args$args[2], quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)
