#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--genes_file", type="character", default=NULL,
              help="File of Ensembl gene ids to subset heatmap to [default %default]" ),
  make_option("--detct", action="store_true", default=FALSE,
              help="Data is DeTCT data, not RNA-seq [default %default]"),
  make_option("--transform", action="store", type="character", default="none",
              help="Transformation to apply to data [default %default]"),
  make_option("--cluster", action="store", type="character", default="none",
              help="Cluster rows, columns, both or none [default %default]"),
  make_option("--colour_palette", type="character", default='plasma',
              help="palette for colour scale [default %default]" ),
  make_option("--cell_colour", type="character", default=NULL,
              help="colour for the outlines of the cells [default %default]" ),
  make_option("--gene_names", action="store_true", type="logical", default=FALSE,
              help="print gene names [default %default]" ),
  make_option("--sample_names", action="store_true", type="logical", default=FALSE,
              help="print sample names [default %default]" ),
  make_option("--metadata_file", type="character", default=NULL,
              help="feather file of metadata to plot with heatmap [default %default]" ),
  make_option("--metadata_ycol", action="store", default='Category',
              help="Column in metadata to plot on y-axis [default %default]"),
  make_option("--metadata_fill", action="store", default='Value',
              help="Column in metadata to plot as fill colour [default %default]"),
  make_option("--metadata_fill_palette", action="store", default=NULL,
              help="Fill palette for metadata plot [default colour-blind friendly palette from biovisr]"),
  make_option("--relative_plot_sizes", action="store", default=NULL,
              help=paste("The relative sizes of the plots as a comma-separated list",
                         "will be treated as heights [default %default]")),
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

# cmd_line_args <- list(
#   options = list(
#     "genes_file" = "no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-cell_pellet/pca/PC1/PC1-pos-genes.tsv",
#     "metadata_file" = "output/sample2organism-paxgene.ftr",
#     "metadata_ycol" = "Organism",
#     "metadata_fill" = "InfectionType",
#     "metadata_fill_palette" = "fill_colour",
#     "detct" = FALSE,
#     "transform" = "center_scale",
#     "cluster" = "both",
#     "colour_palette" = 'plasma',
#     "width" = 7,
#     "height" = 10,
#     "fill_limits" = NULL,
#     "output_count_file" = NULL,
#     "output_rda_file" = NULL
#   ),
#   args = c(
#     'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-paxgene-icu_vs_hv/samples.tsv',
#     'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-paxgene-icu_vs_hv/sig.tsv',
#     'no_icu_004-026_cell-icu_008-062_pax/deseq2-noHb-paxgene-icu_vs_hv/heatmap-genes_clustered-sig.pdf'
#   )
# )

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'gene_expr_heatmap.R',
    usage = "Usage: %prog [options] sample_file count_file output_file" ),
  positional_arguments = 3
)
output_file <- cmd_line_args$args[3]
fill_palette <- cmd_line_args$options[['colour_palette']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
if (!is.null(cmd_line_args$options[['fill_limits']])) {
    fill_limits <-
        strsplit(cmd_line_args$options[['fill_limits']], ",", fixed = TRUE)[[1]]
} else {
    fill_limits <- cmd_line_args$options[['fill_limits']]
}

packages <- c('tidyverse', 'viridis', 'rnaseqtools', 'biovisr', 'patchwork')
if (grepl('svg$', output_file)) { packages <- c(packages, 'svglite') }
if (!is.null(cmd_line_args$options[['metadata_file']])) {
    if (grepl("\\.ftr$", cmd_line_args$options[['metadata_file']])) { packages <- c(packages, 'feather') }
}
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# unpack relative_plot_sizes options
if (is.null(cmd_line_args$options[['relative_plot_sizes']])) {
  relative_plot_sizes <- c(9,1)
} else {
  relative_plot_sizes <- as.integer(unlist(str_split(cmd_line_args$options[['relative_plot_sizes']], ",")))
}

# load sample data
samples <- read_tsv(cmd_line_args$args[1])

#data <- read.delim(cmd_line_args$args[2], check.names = FALSE)
data <- load_rnaseq_data(cmd_line_args$args[2])

# make unique rownames depending on DeTCT/RNA-seq 
if (cmd_line_args$options[['detct']]) {
    data$id <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    data$id <- data[['GeneID']]
}

# subset data to gene list if provided
if (!is.null(cmd_line_args$options[['genes_file']])) {
    genes <- read.delim(file = cmd_line_args$options[['genes_file']], header = FALSE)
    # subset data to genes and arrange it in same order
    data <- filter(data, GeneID %in% genes[[1]]) %>% 
      arrange(match(GeneID, genes[[1]]))
}

# get normalised count data
counts <- get_counts(data, samples, normalised = TRUE)

if (cmd_line_args$options[['transform']] == 'center_scale') {
    counts <- as.data.frame(t( scale( t(counts) ) ))
}

# add gene ids to rownames
rownames(counts) <- data$GeneID

# skip clustering if only 1 row
if (nrow(counts) > 1) {
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
}

# add id column
counts <- counts %>%
  rownames_to_column(., var = "id") %>%
  mutate(id = factor(id, levels = id))

# melt count df
plot_data <- counts %>%
    pivot_longer(-id, names_to = 'Sample', values_to = 'Count')
# set levels of Sample
plot_data$Sample <- factor(plot_data$Sample, levels = unique(plot_data$Sample))

# if (cmd_line_args$options[['sample_names']]) {
#   heatmap_theme <- heatmap_theme + theme(axis.text.x = element_text(colour = "black", angle = 90))
# }

if (is.null(cmd_line_args$options[['cell_colour']]) ) {
  heatmap_plot <- ggplot() +
    geom_tile(data = plot_data,
              aes(x = Sample, y = id, fill = Count))
} else {
  heatmap_plot <- ggplot() +
    geom_tile(data = plot_data, colour = cmd_line_args$options[['cell_colour']],
              aes(x = Sample, y = id, fill = Count))
}

distiller_palettes <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral",
                        "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3",
                        "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", 
                        "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
if (fill_palette %in% c('viridis', 'plasma', 'magma', 'inferno', 'cividis')) {
  heatmap_plot <- heatmap_plot +
    scale_fill_viridis(option = fill_palette, limits = fill_limits)
} else if (fill_palette %in% distiller_palettes) {
  heatmap_plot <- heatmap_plot +
    scale_fill_distiller(palette = fill_palette, limits = fill_limits)
} else {
  stop("fill_palette is not recognised")
}

if (cmd_line_args$options[['gene_names']]) {
  # make function to return a gene name for an Ensmbl id
  ids2names <- function(ids) {
    map_chr(ids, function(id, data){ data$Name[ data$GeneID == id] }, data)
  }
  heatmap_plot <- heatmap_plot + 
    scale_y_discrete(name = NULL, labels = ids2names)
}
heatmap_plot <- heatmap_plot + 
        biovisr::theme_heatmap(xaxis_labels = cmd_line_args$options[['gene_names']], 
                               yaxis_labels = cmd_line_args$options[['sample_names']],
                               base_size = 12) +
        NULL

# load metadata if provided
# The first column should be the sample ids/names (x axis, matching the samples file)
# the y axis and fill variables are specified by metadata_ycol and metadata_fill
if (!is.null(cmd_line_args$options[['metadata_file']])) {
  if (grepl("\\.ftr$", cmd_line_args$options[['metadata_file']])) {
    metadata_for_plot <- read_feather(cmd_line_args$options[['metadata_file']])
  } else {
    metadata_for_plot <- read_tsv(cmd_line_args$options[['metadata_file']])
  }

  # subset metadata to samples
  metadata_for_plot <- filter(metadata_for_plot, sample %in% levels(plot_data$Sample))
  # order sample levels by clustering
  metadata_for_plot[['sample']] <-
    factor(metadata_for_plot[['sample']],
           levels = levels(plot_data$Sample))
  # reverse levels of y axis to make plot match
  category_col <- cmd_line_args$options[['metadata_ycol']]
  if (class(metadata_for_plot[[category_col]]) == 'character') {
    metadata_for_plot[[category_col]] <-
      factor(metadata_for_plot[[category_col]],
              levels = rev(unique(metadata_for_plot[[category_col]])))
  } else if (class(metadata_for_plot[[category_col]]) == 'factor') {
    metadata_for_plot[[category_col]] <-
      factor(metadata_for_plot[[category_col]],
             levels = rev(levels(metadata_for_plot[[category_col]])))
  }
  
  # sort out fill_palette
  fill_col <- cmd_line_args$options[['metadata_fill']]
  # set levels of fill column
  metadata_for_plot[[fill_col]] <-
    factor(metadata_for_plot[[fill_col]],
            levels = unique(metadata_for_plot[[fill_col]]))
  
  fill_palette <- cmd_line_args$options[['metadata_fill_palette']]
  if (!is.null(fill_palette)) {
    if (fill_palette %in% colnames(metadata_for_plot)) {
      # make a named vector of levels to colours
      if( length(unique(metadata_for_plot[[fill_palette]])) ==
                          nlevels(metadata_for_plot[[fill_col]]) ) {
        cat2colour <- metadata_for_plot[ , c(fill_col, fill_palette)] %>% unique()
        fill_palette <- cat2colour[[fill_palette]]
        names(fill_palette) <- cat2colour[[fill_col]]
      }
    }
  }
  
  metadata_plot <-
    df_heatmap(metadata_for_plot, x = 'sample', y = category_col,
               fill = fill_col, fill_palette = fill_palette,
               colour = "grey50", size = 0.5,
               xaxis_labels = FALSE, yaxis_labels = TRUE,
               na.translate = FALSE
    ) + guides(fill = guide_legend(title = fill_col, reverse = FALSE)) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(colour = "black"))
}

if (grepl('ps$', output_file)) {
  postscript(file = output_file, paper = "special", horizontal = FALSE,
              width = plot_width, height = plot_height )
} else if (grepl('svg$', output_file)) {
  svglite(file = output_file, width = plot_width, height = plot_height )
} else { # default to pdf 
  pdf(file = output_file, width = plot_width, height = plot_height )
}
if (is.null(cmd_line_args$options[['metadata_file']])) {
  print(heatmap_plot)
} else {
  print(heatmap_plot + metadata_plot +
          plot_layout(ncol = 1, nrow = 2,
                      heights = relative_plot_sizes))
}
dev.off()

# reverse counts matrix so it fits the heatmap
counts <- counts[ seq(nrow(counts), 1), ]

if (!is.null(cmd_line_args$options[['output_count_file']])) {
    write.table(counts, file = cmd_line_args$options[['output_count_file']],
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

if (!is.null(cmd_line_args$options[['output_rda_file']])) {
  object_list <- c('data', 'counts', 'plot_data', 'heatmap_plot')
  if (!is.null(cmd_line_args$options[['metadata_file']])) {
    object_list <- c(object_list, 'metadata_plot')
  }
  save(list = object_list, file = cmd_line_args$options[['output_rda_file']])
}
