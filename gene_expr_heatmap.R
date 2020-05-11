#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--detct", action="store_true", default=FALSE,
              help="Data is DeTCT data, not RNA-seq [default %default]"),
  make_option("--transform", action="store", type="character", default="none",
              help="Transformation to apply to data [default %default]"),
  make_option("--cluster", action="store", type="character", default="none",
              help="Cluster rows, columns, both or none [default %default]"),
  make_option("--colour_palette", type="character", default='plasma',
              help="palette for colour scale [default %default]" ),
  make_option("--width", type="numeric", default=7,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=10,
              help="height of plot (inches) [default %default]" ),
  make_option("--fill_limits", type="character", default=NULL,
              help="Limits of fill scale (comma-separated list of two numbers) [default %default]" ),
  make_option("--output_count_file", type="character", default=NULL,
              help="Name of file for the transformed and clustered counts [default %default]" ),
  make_option("--output_rda_file", type="character", default=NULL,
              help="output an .rda file containing the plot object and clustered counts [default %default]" )
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'gene_expr_heatmap.R',
    usage = "Usage: %prog [options] sample_file count_file output_file" ),
  positional_arguments = 3
)
output_file <- cmd_line_args$args[3]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
if (!is.null(cmd_line_args$options[['fill_limits']])) {
    fill_limits <-
        strsplit(cmd_line_args$options[['fill_limits']], ",", fixed = TRUE)[[1]]
} else {
    fill_limits <- cmd_line_args$options[['fill_limits']]
}

packages <- c('tidyverse', 'viridis', 'rnaseqtools')
if (grepl('svg$', output_file)) { packages <- c(packages, 'svglite') }
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load sample data
samples <- read.delim(cmd_line_args$args[1], row.names = 1)
# set levels of condition to order in which they appear in the file
samples$condition <- factor(samples$condition,
                            levels = unique(samples$condition))

#data <- read.delim(cmd_line_args$args[2], check.names = FALSE)
data <- load_rnaseq_data(cmd_line_args$args[2])

## Support different column names
#names(data)[names(data) == 'chr']               <- 'Chr'
#names(data)[names(data) == '#Chr']              <- 'Chr'
#names(data)[names(data) == 'start']             <- 'Start'
#names(data)[names(data) == 'end']               <- 'End'
#names(data)[names(data) == 'strand']            <- 'Strand'
#names(data)[names(data) == 'ID']                <- 'Gene ID'
#names(data)[names(data) == 'adjpval']           <- 'adjp'
#names(data)[names(data) == 'padj']              <- 'adjp'
#names(data)[names(data) == 'Adjusted p value']  <- 'adjp'
#names(data)[names(data) == 'Gene name']         <- 'Name'
#names(data)[ grepl("e[0-9]+ Ensembl Gene ID", names(data)) ] <- 'Gene ID'
#names(data)[names(data) == 'GeneID']              <- 'Gene ID'

# 
if (cmd_line_args$options[['detct']]) {
    rownames(data) <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    rownames(data) <- data[['GeneID']]
}

counts <- data[ , grepl("normalised.count$", names(data)) ]
names(counts) <- sub(".normalised.count", "", names(counts))

# make sure counts samples are in same order as samples
counts <- counts[ , rownames(samples) ]

if (cmd_line_args$options[['transform']] == 'center_scale') {
    counts <- as.data.frame(t( scale( t(counts) ) ))
}

# Clustering
cluster <- function(df) {
    df <- as.data.frame(df)
    distance_matrix <- dist(df)
    clustering <- hclust(distance_matrix)
    df_ordered <- df[ clustering$order, ]
    return(df_ordered)
}

cluster_rows <- FALSE
cluster_columns <- FALSE
if (cmd_line_args$options[['cluster']] == "rows") {
    cluster_rows <- TRUE
} else if (cmd_line_args$options[['cluster']] == "columns") {
    cluster_columns <- TRUE
} else if (cmd_line_args$options[['cluster']] == "both") {
    cluster_rows <- TRUE
    cluster_columns <- TRUE
}
if (cluster_rows) {
    counts <- cluster(counts)
}
if (cluster_columns) {
    counts <- as.data.frame(t( cluster(t(counts)) ))
}

counts$id <- factor(rownames(counts), levels = rownames(counts))

# melt count df
plot_data <- counts %>%
    pivot_longer(-id, names_to = 'Sample', values_to = 'Count')
# set levels of Sample
plot_data$Sample <- factor(plot_data$Sample, levels = unique(plot_data$Sample))

heatmap_plot <- ggplot() +
        geom_raster(data = plot_data, aes(x = Sample, y = id,
                                           fill = Count)) +
        scale_fill_viridis(option = cmd_line_args$options[['colour_palette']],
                            limits = fill_limits) +
        theme_void() +
        NULL

if (grepl('ps$', output_file)) {
  postscript(file = output_file, paper = "special", horizontal = FALSE,
              width = plot_width, height = plot_height )
} else if (grepl('svg$', output_file)) {
  svglite(file = output_file, width = plot_width, height = plot_height )
} else { # default to pdf 
  pdf(file = output_file, width = plot_width, height = plot_height )
}
print(heatmap_plot)
dev.off()

# reverse counts matrix so it fits the heatmap
counts <- counts[ seq(nrow(counts), 1),
                c('id', colnames(counts)[ !grepl("id", colnames(counts))] )]

if (!is.null(cmd_line_args$options[['output_count_file']])) {
    write.table(counts, file = cmd_line_args$options[['output_count_file']],
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

if (!is.null(cmd_line_args$options[['output_rda_file']])) {
    save(counts, heatmap_plot, file = cmd_line_args$options[['output_rda_file']])
}
