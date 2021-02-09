# Script to plot a volcano plot from DESeq2 results

# load packages
library('optparse')

option_list <- list(
  make_option("--labels", action="store_true", type="logical", default=FALSE,
              help="Label the points which are above the logfc and pvalue thresholds [default %default]" ),
  make_option("--log2fc_threshold", type="numeric", default=NULL,
              help="Log2 FC threshold above which to label the points [default %default]" ),
  make_option("--pval_threshold", type="numeric", default=NULL,
              help="P value cut-off to use to colour/label the points [default %default]" )
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'volcano_plot.R',
    usage = "Usage: %prog [options] results_file output_file",
    description = 'This script makes a volcano plot from a DESeq2 results file'),
  positional_arguments = 2
)

# cmd_line_args <-
#    list(options = list(),
#        args = c('deseq2/3dpf_uninf_hom_vs_3dpf_uninf_sib.tsv',
#                 'deseq2/3dpf_uninf_hom_vs_3dpf_uninf_sib.volcano.pdf'))

# unpack options
deseq_results_file <- cmd_line_args$args[1]
output_file <- cmd_line_args$args[2]

# load packages 
packages <- c('tidyverse', 'viridis', 'ggrepel', 'ggrastr')
if (grepl('svg$', output_file)) { packages <- c(packages, 'svglite') }
for( package in packages ){
    suppressPackageStartupMessages(
        suppressWarnings( library(package, character.only = TRUE) ) )
}

# load data
deseq_results <- read.delim(deseq_results_file) %>% 
  filter(., !is.na(adjp)) # remove gene with adjp == NA

# make -log10p column
deseq_results$log10p <- -log10(deseq_results$adjp)

# make factor for up and down genes
deseq_results$sig <- deseq_results$adjp < 0.05
deseq_results$up_or_down <- rep(NA, nrow(deseq_results))
deseq_results$up_or_down[ deseq_results$sig & deseq_results$log2fc > 0 ] <- 'up'
deseq_results$up_or_down[ deseq_results$sig & deseq_results$log2fc < 0 ] <- 'down'
deseq_results$up_or_down <- factor(deseq_results$up_or_down,
                                   levels = c('up', 'down'))

# sort by Adjusted pvalue
deseq_results <- arrange(deseq_results, desc(adjp))

# highlight specific genes
# label genes with abs(fold change) > 2
if (cmd_line_args$options[['labels']]) {
  # set defaults if not set
  if (is.null(cmd_line_args$options[['log2fc_threshold']]) &
      is.null(cmd_line_args$options[['pval_threshold']])) {
    cmd_line_args$options[['log2fc_threshold']] <- 2
    cmd_line_args$options[['pval_threshold']] <- 1e-5
  }
  if (!is.null(cmd_line_args$options[['log2fc_threshold']])) {
    biggest_changers <- 
      filter(deseq_results, 
             abs(log2fc) > cmd_line_args$options[['log2fc_threshold']])
  } else {
    biggest_changers <- deseq_results
  }
  if (!is.null(cmd_line_args$options[['pval_threshold']])) {
    biggest_changers <- 
      filter(deseq_results, 
             adjp < cmd_line_args$options[['pval_threshold']] )
  }
}

# plot adjusted pvalue against log2 fold change
# with up coloured in orange and down coloured in blue
# label highest changers
volcano_plot <-
    ggplot(data = deseq_results, aes(x = log2fc, y = log10p, colour = up_or_down)) +
        rasterise(geom_point(), dpi = 300)
if (cmd_line_args$options[['labels']]) {
  volcano_plot <- volcano_plot +
        geom_text_repel(data = biggest_changers, aes(label = Name)) 
}

volcano_plot <- volcano_plot +
  scale_colour_manual(values = c(up = '#cc6600', down = '#0073b3'), na.value = "grey80") +
  theme_minimal() +
  guides(colour = "none") +
  labs(x = expr(log[2]*'(Fold Change)'), y = expr(log[10]*'(Adjusted pvalue)')) +
  NULL

# output plot
# get output type from filename suffix
# pdf is default if nothing matches
if (sub("^.*\\.", "", output_file) == "eps") {
  postscript(file = output_file)
} else if (sub("^.*\\.", "", output_file) == "svg") {
  library('svglite')
  svglite(file = output_file)
} else {
  pdf(file = output_file)
}

# print plot and close graphics device
print(volcano_plot)
dev.off()
