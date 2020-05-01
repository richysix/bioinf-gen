#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-o", "--output_file"), action="store", default='plots/sample_clustering.pdf',
              help="Name of output file [default]"),
  make_option("--plots_rda_file", action="store", default=NULL,
              help="Name of output file [default]"),
  make_option("--metadata_file", action="store", default=NULL,
              help="Name of output file [default]"),
  make_option(c("-d", "--debug"), action="store_true", default=FALSE,
              help="Print debugging info [default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

desc <- paste('\nScript to cluster samples from an RNA-Seq experiment', sep = "\n")

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'cluster_samples_rnaseq.R',
    description = desc,
    usage = "Usage: %prog [options] all_file samples_file" ),
  positional_arguments = 2
)

# cmd_line_args <- list(
#   options = list("output_file" = 'test-cor-plot.png',
#                  "metadata_file" = 'test-metadata.tsv',
#                  "debug" = TRUE, "verbose" = TRUE),
#   args = c('all.tsv', 'test-samples.tsv')
# )

packages <- c('rnaseqtools', 'biovisr', 'miscr', 'tidyverse', 'patchwork')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load RNAseq data
data <- load_rnaseq_data(datafile = cmd_line_args$args[1])

# load sample data
samples <- read_tsv(cmd_line_args$args[2])

# get count data
get_counts <- function(data, samples = NULL) {
  count_data <- data[,grepl(".normalised.counts?$", names(data))]
  names(count_data) <- gsub(".normalised.counts?$", "", names(count_data))
  
  # Subset and reorder count data
  # TO ADD: check all samples exist in counts
  # and if any are in counts that aren't in samples
  if (!is.null(samples)) {
    count_data <- count_data[, samples$sample ]
  }
  
  return(count_data)
}
count_data <- get_counts(data, samples)

# cluster samples
# POSSIBLE OPTION: method	(the agglomeration method to use)
# the cluster function will pass that through to hclust
clustering <- cluster(as.matrix(count_data), scale = FALSE, 
                      dist_method = "spearman", clustering = TRUE)

# plot tree
tree_plot <- dendro_plot(clustering$clustering)

# load metadata if it exists
# The first column should be the sample ids/names (x axis, matching the samples file)
# the second column are the metadata catgories (y axis)
# the third column are the values (fill)
if (!is.null(cmd_line_args$options[['metadata_file']])) {
  metadata <- read_tsv(cmd_line_args$options[['metadata_file']])
} else {
  metadata <- NULL
}
# metadata plot
if (!.is.null(metadata)) {
  # set order of samples
  y_cat <- names(metadata)[1]
  metadata <- select(metadata, y_cat, colnames(clustering$matrix)) %>% 
    pivot_longer(., -!!y_cat, names_to = "sample", values_to = "Category")
  metadata[[y_cat]] <- factor(metadata[[y_cat]],
                              levels = rev(unique(metadata[[y_cat]])))
  
  metadata[['sample']] <- factor(metadata[['sample']],
                                 levels = unique(metadata[['sample']]))
  
  metadata[['Category']] <- factor(metadata[['Category']],
                                   levels = rev(unique(metadata[['Category']])))
  
  metadata_plot <- 
    df_heatmap(metadata, x = 'sample', y = y_cat, 
               fill = 'Category', colour = "black", size = 0.8,
               xaxis_labels = FALSE, yaxis_labels = FALSE,
               na.translate = FALSE, fill_palette = NULL
    ) + theme(axis.title.y = element_blank())
}

# plot cor matrix
# create correlation matrix and reorder
# as clustering
cor_matrix <- cor(as.matrix(count_data))
cor_matrix <- cor_matrix[ , clustering$clustering$order ]

cor_matrix_plot <- 
  matrix_heatmap(cor_matrix, x_title = "Sample", y_title = "Sample",
                 fill_title = "Spearman", fill_palette = "plasma",
                 xaxis_labels = TRUE, yaxis_labels = TRUE,
                 base_size = 14)

# output plots together
postscript(file = cmd_line_args$options[['output_file']], 
            width = 7, height = 10)
if (is.null(metadata)) {
  print(tree_plot + cor_matrix_plot + plot_layout(ncol = 1, nrow = 2))
} else {
  print(tree_plot + metadata_plot + cor_matrix_plot + 
          plot_layout(ncol = 1, nrow = 3, heights = c(4,1,5)))
}
dev.off()

if (!is.null(cmd_line_args$options[['plots_rda_file']])) {
  save(list(tree_plot, metadata_plot, cor_matrix_plot), 
       file = cmd_line_args$options[['plots_rda_file']])
}

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2020. University of Cambridge.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
  