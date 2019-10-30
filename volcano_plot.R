# Script to plot a volcano plot from DESeq2 results

# load packages
library('optparse')

option_list <- list(
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'volcano_plot.R',
    usage = "Usage: %prog [options] results_file output_file",
    description = 'This script makes a volcano plot from a DESeq2 results file'),
  positional_arguments = 3
)

#cmd_line_args <-
#    list(options = list(),
#        args = c('deseq2/3dpf_uninf_hom_vs_3dpf_uninf_sib.tsv',
#                 'deseq2/3dpf_uninf_hom_vs_3dpf_uninf_sib.volcano.pdf'))

# unpack options
deseq_results_file <- cmd_line_args$args[1]
output_file <- cmd_line_args$args[2]

# load packages 
packages <- c('tidyverse', 'viridis', 'ggrepel')
if (grepl('svg$', output_file)) { packages <- c(packages, 'svglite') }
for( package in packages ){
    suppressPackageStartupMessages(
        suppressWarnings( library(package, character.only = TRUE) ) )
}

# load data
deseq_results <- read.delim(deseq_results_file)

# make -log10p column
deseq_results$log10p <- -log10(deseq_results$adjp)

# make factor for up and down genes
deseq_results$sig <- !is.na(deseq_results$adjp) & deseq_results$adjp < 0.05
deseq_results$up_or_down <- rep(NA, nrow(deseq_results))
deseq_results$up_or_down[ deseq_results$sig & deseq_results$log2fc > 0 ] <- 'up'
deseq_results$up_or_down[ deseq_results$sig & deseq_results$log2fc < 0 ] <- 'down'
deseq_results$up_or_down <- factor(deseq_results$up_or_down,
                                   levels = c('up', 'down'))

# sort by Adjusted pvalue
deseq_results <- arrange(deseq_results, desc(adjp))

# highlight specific genes
# label genes with abs(fold change) > 2
biggest_changers <- filter(deseq_results, abs(log2fc) > 3)

# plot adjusted pvalue against log2 fold change
# with up coloured in orange and down coloured in blue
# label highest changers
volcano_plot <-
    ggplot(data = deseq_results, aes(x = log2fc, y = log10p, colour = up_or_down)) +
        geom_point() +
        geom_text_repel(data = biggest_changers, aes(label = Name)) +
        scale_colour_manual(values = c(up = '#cc6600', down = '#0073b3'), na.value = "grey80") +
        theme_minimal() +
        guides(colour = "none") +
        labs(x = expr(log[2]*'(Fold Change)'), y = expr(log[10]*'(Adjusted pvalue)')) +
        NULL


