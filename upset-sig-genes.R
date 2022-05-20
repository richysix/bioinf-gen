#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--plot_file", type="character", default='upset-plot.tsv',
              help="Name for UpSet plot file [default %default]" ),
  make_option("--output_file", type="character", default='upset-sets.tsv',
              help="Name for output file [default %default]" ),
  make_option("--keep_order", action="store_true", type="logical", default=FALSE,
              help="Keep the order of the sets as on the command line [default %default]" ),
  make_option("--width", type="numeric", default=10,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=8,
              help="height of plot (inches) [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to find and plot the intersections between DE gene lists',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'upset-sig-genes.R',
    usage = "Usage: %prog [options] input_files",
    description = desc ),
  positional_arguments = c(2, Inf)
)

# load packages
packages <- c('tidyverse', 'UpSetR', 'rlang', 'miscr', 'rnaseqtools')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# unpack options
output_file <- cmd_line_args$options[['output_file']]
plot_file <- cmd_line_args$options[['plot_file']]
keep_order <- cmd_line_args$options[['keep_order']]
width <- cmd_line_args$options[['width']]
height <- cmd_line_args$options[['height']]
debug <- cmd_line_args$options[['debug']]

# open input files
sig_genes_list <- purrr::map(cmd_line_args$args, load_rnaseq_data)
# label sets by removing deseq2- from the start and [/.]sig.tsv from the end of the filename
sets <- sub("^deseq2-", "", cmd_line_args$args) %>% 
  sub("[/.]sig.tsv$", "", .)
names(sets) <- sets
names(sig_genes_list) <- sets

# create binary matrix of gene ids to Sets
# make data frame of gene ids and the sets they are present in 
gene_ids_binary_matrix <- purrr::map_dfr(sig_genes_list, function(x){
  select(x, GeneID) %>% 
    mutate(., Present = 1)
  },
  .id = "Set") %>% 
  # then turn that into a binary matrix with columns
  # for the sets and 1/0 to indicate presence/absence in that set
  pivot_wider(., id_cols = GeneID, names_from = Set,
              values_from = Present, values_fill = 0)

# Rearrange sets by size unless --keep_order is true
if (!keep_order) {
  set_order <- purrr::map_dbl(sets, ~ sum(gene_ids_binary_matrix[[.x]])) %>% 
    sort(., decreasing = TRUE) %>% names()
  gene_ids_binary_matrix <- select(gene_ids_binary_matrix, GeneID, all_of(set_order))
}

# work out which intersection each gene is in
# each intersection is represented as a number whose binary representation indicates which set it is a part of
# e.g. for 3 sets, 5 = 101 is a gene in Set3 and Set1
gene_ids_to_sets <- gene_ids_binary_matrix %>% 
  unite("Intersection", !matches("GeneID"), remove = FALSE, sep = "") %>% 
  mutate(Intersection = strtoi(Intersection, base = 2))

# convert the 1/0 columns to strings
convert_to_string <- function(df, col_name) {
  col_name_sym <- rlang::sym(col_name)
  df %>% mutate(!!col_name_sym := case_when(
    !!col_name_sym == 0 ~ NA_character_,
    !!col_name_sym == 1 ~ col_name
  ))
}
for (col_name in sets) {
  gene_ids_to_sets <- convert_to_string(gene_ids_to_sets, col_name)
}
# stitch the strings together to a string representation of each intersection
# and join in Gene Names
gene_ids_to_sets <- gene_ids_to_sets %>% 
  unite("Set", all_of(sets), sep = "|", na.rm = TRUE) %>% 
  mutate(., Intersection = factor(Intersection),
         Intersection = fct_infreq(Intersection)) %>% # order by the size of the intersections
  arrange(., Intersection)

# output to file
write_tsv(gene_ids_to_sets, file = output_file)

# list of genes present in each set
gene_ids_list <- purrr::map(sig_genes_list, ~ .x$GeneID)

# do upset diagram
open_graphics_device(filename = plot_file,
                     width = width, height = height)
upset(fromList(gene_ids_list), order.by = "freq",
      sets = sets, keep.order = keep_order)
invisible(dev.off())

# options for testing
# cmd_line_args <-
#   list(options = list('output_file' = 'upset-test.tsv',
#                       'plot_file' = 'upset-test-2.svg',
#                       'keep_order' = FALSE,
#                       'width' = 10, 'height' = 8,
#                       debug = TRUE),
#        args = c('deseq2-srpk3_hom_ttnb_het_vs_srpk3_hom_ttnb_wt/sig.tsv',
#                 'deseq2-srpk3_hom_ttnb_het_vs_srpk3_wt_ttnb_het/sig.tsv',
#                 'deseq2-srpk3_hom_ttnb_het_vs_srpk3_wt_ttnb_wt/sig.tsv',
#                 'deseq2-srpk3_hom_ttnb_wt_vs_srpk3_wt_ttnb_wt/sig.tsv',
#                 'deseq2-srpk3_wt_ttnb_het_vs_srpk3_wt_ttnb_wt/sig.tsv'))
