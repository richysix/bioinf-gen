#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-o", "--output_file"), action="store", default='sample_clustering.svg',
              help="Name of output file [default %default]"),
  make_option("--plots_rda_file", action="store", default=NULL,
              help="Name of output file [default %default]"),
  make_option("--distance_measure", action="store", default='spearman',
              help="Measure to use to calculate the distance matrix [default %default]"),
  make_option("--metadata_file", action="store", default=NULL,
              help="Name of output file [default %default]"),
  make_option("--plot_width", action="store", default=12,
              help="Width of plot [default %default]"),
  make_option("--plot_height", action="store", default=20,
              help="Height of plot [default %default]"),
  make_option(c("-d", "--debug"), action="store_true", default=FALSE,
              help="Print debugging info [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default %default]")
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
#   options = list("output_file" = 'no_outliers-icu_004_cell_pellet-icu_026_cell_pellet/deseq2-cell_pellet/sample_cor.svg',
#                  "plots_rda_file" = 'no_outliers-icu_004_cell_pellet-icu_026_cell_pellet/deseq2-cell_pellet/sample_cor.rda',
#                  "metadata_file" = 'output/infection-status-cell_pellet.ftr',
#                  "debug" = TRUE, "verbose" = TRUE),
#   args = c('no_outliers-icu_004_cell_pellet-icu_026_cell_pellet/deseq2-cell_pellet/all.tsv',
#            'no_outliers-icu_004_cell_pellet-icu_026_cell_pellet/deseq2-cell_pellet/samples.tsv')
# )

packages <- c('rnaseqtools', 'biovisr', 'miscr', 'tidyverse', 'patchwork',
              'feather', 'svglite')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# set verbose options
if (!cmd_line_args$options[['verbose']]){
  options(readr.show_progress=FALSE)
}

# load RNAseq data
data <- load_rnaseq_data(data_file = cmd_line_args$args[1])

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
                      dist_method = cmd_line_args$options[['distance_measure']], 
                      clustering = TRUE)

# plot tree
tree_plot <- dendro_plot(clustering$clustering)

# load metadata if it exists
# The first column should be the sample ids/names (x axis, matching the samples file)
# the second column are the metadata catgories (y axis)
# the third column are the values (fill)
if (!is.null(cmd_line_args$options[['metadata_file']])) {
  if (grepl("\\.ftr$", cmd_line_args$options[['metadata_file']])) {
    metadata_for_plot <- read_feather(cmd_line_args$options[['metadata_file']])
  } else {
    metadata_for_plot <- read_tsv(cmd_line_args$options[['metadata_file']])
  }
} else {
  metadata_for_plot <- NULL
}
# metadata plot
if (!is.null(metadata_for_plot)) {
  # subset data to samples in samples file
  metadata_for_plot <- filter(metadata_for_plot, sample %in% samples$sample)
  # order levels by clustering
  metadata_for_plot[['sample']] <- 
    factor(metadata_for_plot[['sample']],
           levels = colnames(clustering$matrix))
  # reverse levels of Category to make plot match
  metadata_for_plot[['Category']] <- factor(metadata_for_plot[['Category']],
                                   levels = rev(levels(metadata_for_plot[['Category']])))

  metadata_plot <- 
    df_heatmap(metadata_for_plot, x = 'sample', y = "Category", 
               fill = "Value", colour = "black", size = 0.8,
               xaxis_labels = FALSE, yaxis_labels = FALSE,
               na.translate = FALSE, fill_palette = NULL
    ) + theme(axis.title.y = element_blank(),
              axis.text.y = element_text(colour = "black"))
}

# plot cor matrix
# create correlation matrix and reorder
# as clustering
cor_matrix <- cor(as.matrix(count_data))
cor_matrix <- cor_matrix[ clustering$clustering$order , clustering$clustering$order ]

cor_matrix_plot <- 
  matrix_heatmap(cor_matrix, x_title = "Sample", y_title = "Sample",
                 fill_title = "Spearman", fill_palette = "plasma",
                 xaxis_labels = TRUE, yaxis_labels = TRUE,
                 base_size = 10)

# output plots together
svglite(file = cmd_line_args$options[['output_file']], 
        width = cmd_line_args$options[['plot_width']], 
        height = cmd_line_args$options[['plot_height']])
if (is.null(metadata_for_plot)) {
  print(tree_plot + cor_matrix_plot + plot_layout(ncol = 1, nrow = 2))
} else {
  print(tree_plot + metadata_plot + cor_matrix_plot + 
          plot_layout(ncol = 1, nrow = 3, heights = c(4,1,5)))
}
dev.off()

if (!is.null(cmd_line_args$options[['plots_rda_file']])) {
  save(tree_plot, metadata_plot, cor_matrix_plot, 
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
  