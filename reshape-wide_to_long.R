#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--debug", type="logical", default=FALSE, action="store_true",
            help="Turns on debugging statements [default %default]" ),
  make_option("--columns", type="character", default=NULL,
              help="Name of variables to pivot [default All but the first column]" ),
  make_option("--variable.name", type="character", default='variable',
              help="Name of variable for measured variable names [default %default]" ),
  make_option("--value.name", type="character", default='value',
              help="Name of variable used to store values [default %default]" )
)

desc <- paste('', 'Script to reformat a file from wide to long data',
    sep = "\n")

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'reshape-wide_to_long.R',
    usage = "Usage: %prog [options] inputFile", description = desc ),
  positional_arguments = 2
)

library('tidyverse')

input_file <- cmd_line_args$args[1]
if (input_file == '-') {
    input_file <- file("stdin")
}
data <- read_tsv(input_file)

# ID variables
if (is.null(cmd_line_args$options[['columns']])) {
  col_ids <- setdiff(seq_len(ncol(data)), c(1))
} else {
  col_names <- unlist(strsplit(cmd_line_args$options[['columns']], ','))
  # check variable names
  col_names_exist <- col_names %in% names(data)
  if (any(!col_names_exist)) {
    missing_cols <- col_names[ !col_names_exist ]
    err_msg <- paste0('At least one of the --columns names are missing from the data\n',
                      'Missing: ', paste0(missing_cols, collapse = ', '), '\n')
    stop(err_msg)
  }
  col_ids <- which(colnames(data) %in% col_names)
}

data <- data %>%
  pivot_longer(., cols = all_of(col_ids), 
               names_to = cmd_line_args$options[["variable.name"]], 
               values_to = cmd_line_args$options[["value.name"]])
  # mutate(., category = factor(category, levels = unique(category))) %>% 
  # arrange(., category, value) %>% 
  # mutate(., value = factor(value, levels = unique(value)))

output_file <- cmd_line_args$args[2]
# if (output_file == '-') {
#   output_file <- file("stdout")
# }
write_tsv(data, output_file)
