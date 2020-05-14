#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-o", "--output_file"), action="store", default='sample_clustering.svg',
              help="Name of output file [default %default]"),
  make_option("--plots_rda_file", action="store", default=NULL,
              help="Name of plots rda file [default %default]"),
  make_option("--distance_measure", action="store", default='spearman',
              help="Measure to use to calculate the distance matrix [default %default]"),
  make_option("--metadata_file", action="store", default=NULL,
              help="Name of metadata file [default %default]"),
  make_option("--metadata_ycol", action="store", default='Category',
              help="Column in metadata to plot on y-axis [default %default]"),
  make_option("--metadata_fill", action="store", default='Value',
              help="Column in metadata to plot as fill colour [default %default]"),
  make_option("--metadata_fill_palette", action="store", default=NULL,
              help="Fill palette for metadata plot [default colour-blind friendly palette from biovisr]"),
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
#   options = list("output_file" = 'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-cell_pellet/sample_cor.svg',
#                  "plots_rda_file" = 'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-cell_pellet/sample_cor.rda',
#                  "metadata_file" = 'output/sample2organism-cell_pellet.ftr',
#                  "metadata_ycol" = 'Organism',
#                  "metadata_fill" = 'InfectionType',
#                  "metadata_fill_palette" = NULL,
#                  "plot_width" = 12, "plot_height" = 30,
#                  "distance_measure" = 'spearman',
#                  "debug" = TRUE, "verbose" = TRUE),
#   args = c('deseq2-noHb/all.tsv',
#            'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-cell_pellet/samples.tsv')
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
# the y axis and fill variables are specified by metadata_ycol and metadata_fill
if (!is.null(cmd_line_args$options[['metadata_file']])) {
  if (grepl("\\.ftr$", cmd_line_args$options[['metadata_file']])) {
    metadata_for_plot <- read_feather(cmd_line_args$options[['metadata_file']])
  } else {
    metadata_for_plot <- read_tsv(cmd_line_args$options[['metadata_file']])
  }
}
# metadata plot
# subset data to samples in samples file
metadata_for_plot <- filter(metadata_for_plot, sample %in% samples$sample)
# order levels by clustering
metadata_for_plot[['sample']] <- 
  factor(metadata_for_plot[['sample']],
         levels = colnames(clustering$matrix))
# reverse levels of y axis to make plot match
ycol <- cmd_line_args$options[['metadata_ycol']]
if (class(metadata_for_plot[[ycol]]) == 'character') {
  metadata_for_plot[[ycol]] <- 
    factor(metadata_for_plot[[ycol]],
            levels = rev(unique(metadata_for_plot[[ycol]])))
} else if (class(metadata_for_plot[[ycol]]) == 'factor') {
  metadata_for_plot[[ycol]] <- 
    factor(metadata_for_plot[[ycol]],
           levels = rev(levels(metadata_for_plot[[ycol]])))
}

# sort out fill_palette
fill_col <- cmd_line_args$options[['metadata_fill']]

fill_palette <- cmd_line_args$options[['metadata_fill_palette']]
if (!is.null(fill_palette)) {
  if (fill_palette %in% colnames(metadata_for_plot)) {
    # make a named vector of levels to colours
    if( length(unique(metadata_for_plot[[fill_palette]])) == 
                        nlevels(metadata_for_plot[[fill_col]]) ) {
      fill_palette <- unique(metadata_for_plot[[fill_palette]])
      names(fill_palette) <- levels(metadata_for_plot[[fill_col]])
    }
  }
}

metadata_plot <- 
  df_heatmap(metadata_for_plot, x = 'sample', y = ycol, 
             fill = fill_col, fill_palette = fill_palette,
             colour = "black", size = 0.8,
             xaxis_labels = FALSE, yaxis_labels = TRUE,
             na.translate = FALSE
  ) + guides(fill = guide_legend(title = fill_col, reverse = FALSE)) +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(colour = "black"))

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
          plot_layout(ncol = 1, nrow = 3, heights = c(3,2,5)))
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
  