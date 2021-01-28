#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--debug", type="logical", default=FALSE, action="store_true",
            help="Turns on debugging statements [default %default]" )
)

desc <- paste('', 'Script to reformat a file from wide to long data',
    'Currently assumes that everything but the sample column is to be pivoted.',
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
metadata <- read_tsv(input_file)

metadata <- metadata %>%
  pivot_longer(., cols = -sample, names_to = "category", values_to = "value") %>% 
  mutate(., category = factor(category, levels = unique(category))) %>% 
  arrange(., category, value) %>% 
  mutate(., value = factor(value, levels = unique(value)))

write_tsv(metadata, cmd_line_args$args[2])
