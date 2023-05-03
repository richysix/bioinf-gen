#!/usr/bin/env Rscript

desc <- "
This script plots a volcano plot (Log2[Fold Change] vs -log10[adjusted p]) from a DESeq2 results file.
It will optionally label points that are above pvalue and/or log2[Fold Change] thresholds.
It expects columns named 'adjp' and 'log2fc' in the input file. 
A column called 'Name' is required for labels.
Available output formats are 'eps', 'svg' and 'pdf'
If the ggrastr package is installed, it is used to rasterise the points layer of the plot to reduce the size of the output file.
"

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2022. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007

# load packages
library('optparse')

option_list <- list(
  make_option("--labels", action="store_true", type="logical", default=FALSE,
              help="Label the points which are above the logfc and pvalue thresholds [default %default]" ),
  make_option("--log2fc_threshold", type="numeric", default=NULL,
              help="Absolute Log2 FC threshold above which to label the points [default %default]" ),
  make_option("--pval_threshold", type="numeric", default=NULL,
              help="P value cut-off to use to colour/label the points [default %default]" ),
  make_option("--log2fc_or_pval", action="store_true", type="logical", default=FALSE,
              help=paste0("Change labelling to be points over either threshold [default %default]\n",
                          "The default behaviour if both thresholds are set is to labels points over both thresholds\n") ),
  make_option("--width", type="numeric", default=NULL,
              help="Width of plot in inches [default %default]" ),
  make_option("--height", type="numeric", default=NULL,
              help="Height of plot in inches [default %default]" ),
  make_option("--text_base_size", type="numeric", default=16,
              help="Base size of text for ggplot theme [default %default]" )
)

# For testing. If running this script interactively the options
# get set to defaults and positional arguments are set to
# whatever is in the arguments vector below
if (any(commandArgs() == "--interactive")) {
  arguments <- c('test_data/volcano-test-data.tsv',
                 'test-volcano.pdf')
} else {
  arguments <- commandArgs(trailingOnly = TRUE)
}

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'volcano_plot.R',
    usage = "Usage: %prog [options] results_file output_file",
    description = desc ),
  args = arguments,
  positional_arguments = 2
)

# unpack options
deseq_results_file <- cmd_line_args$args[1]
output_file <- cmd_line_args$args[2]

# load packages 
packages <- c('tidyverse', 'ggrepel', 'miscr')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) ) )
}

if (!is.null(cmd_line_args$options[['log2fc_threshold']]) |
    !is.null(cmd_line_args$options[['pval_threshold']])) {
  cmd_line_args$options[['labels']] <- TRUE
}

# load data
deseq_results <- read_tsv(deseq_results_file) %>% 
  filter(., !is.na(adjp)) %>% # remove gene with adjp == NA
  mutate(., log10p = -log10(adjp),
         up_or_down = factor(case_when( # make factor for up and down genes
           adjp < 0.05 & log2fc > 0 ~ 'up',
           adjp < 0.05 & log2fc < 0 ~ 'down'
         ), levels = c('up', 'down'))
  ) %>% 
  arrange(., desc(adjp)) # sort by Adjusted pvalue

# highlight specific genes
if (cmd_line_args$options[['labels']]) {
  # set defaults if not set
  if (is.null(cmd_line_args$options[['log2fc_threshold']]) &
      is.null(cmd_line_args$options[['pval_threshold']])) {
    # abs(fold change) > 2
    cmd_line_args$options[['log2fc_threshold']] <- 2
    # adjp < 1e-5
    cmd_line_args$options[['pval_threshold']] <- 1e-5
  }
  deseq_results$name_label <- ""
  if (!is.null(cmd_line_args$options[['log2fc_threshold']])) {
    if (!is.null(cmd_line_args$options[['pval_threshold']])) {
      if (cmd_line_args$options[['log2fc_or_pval']]) {
        to_label <- abs(deseq_results$log2fc) > cmd_line_args$options[['log2fc_threshold']] |
          deseq_results$adjp < cmd_line_args$options[['pval_threshold']]
      } else {
        to_label <- abs(deseq_results$log2fc) > cmd_line_args$options[['log2fc_threshold']] &
          deseq_results$adjp < cmd_line_args$options[['pval_threshold']]
      }
    } else {
      to_label <- abs(deseq_results$log2fc) > cmd_line_args$options[['log2fc_threshold']]
    }
  } else {
    to_label <- deseq_results$adjp < cmd_line_args$options[['pval_threshold']]
  }
  deseq_results$name_label[ to_label ] <- deseq_results$Name[ to_label ]
}

# plot adjusted pvalue against log2 fold change
# with up coloured in orange and down coloured in blue
# label genes above the pvalue/log2fc thresholds
volcano_plot <-
  ggplot(data = deseq_results, aes(x = log2fc, y = log10p, colour = up_or_down))

# Use ggrastr if it is installed to make final image smaller
if(suppressPackageStartupMessages(suppressWarnings(require('ggrastr')))) {
  volcano_plot <- volcano_plot +
    rasterise(geom_point(), dpi = 300)
} else { # otherwise use geom_point
  cat("Package ggrastr not available. Continuing without it.\n")
  volcano_plot <- volcano_plot +
    geom_point()
}

# Add gene name labels
if (cmd_line_args$options[['labels']]) {
  volcano_plot <- volcano_plot +
    geom_text_repel(aes(label = name_label),
                    size = cmd_line_args$options[['text_base_size']]*0.75/.pt,
                    point.padding = 0.25, box.padding = 1,
                    force = 2)
}

# Add colour scale, change theme and add axis labels
volcano_plot <- volcano_plot +
  scale_colour_manual(values = c(up = '#cc6600', down = '#0073b3'), na.value = "grey80") +
  theme_minimal(base_size = cmd_line_args$options[['text_base_size']]) +
  guides(colour = "none") +
  labs(x = expr(log[2]*'(Fold Change)'), y = expr(log[10]*'(Adjusted pvalue)')) +
  NULL

# output plot
width <- cmd_line_args$options[['width']]
height <- cmd_line_args$options[['height']]
if (grepl("png$", output_file)) {
  width <- ifelse(is.null(width), 480, width)
  height <- ifelse(is.null(height), 480, height)
} else {
  width <- ifelse(is.null(width), 7, width)
  height <- ifelse(is.null(height), 7, height)
}

miscr::open_graphics_device(
  filename = output_file,
  width = width,
  height = height
)
print(volcano_plot)
dev.off()
