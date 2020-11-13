#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--genes_file", type="character", default=NULL,
              help="File of Ensembl gene ids to subset heatmap to [default %default]" ),
  make_option("--detct", action="store_true", default=FALSE,
              help="Data is DeTCT data, not RNA-seq [default %default]"),
  make_option("--transform", action="store", type="character", default="rlog",
              help="Transformation to apply to data [default %default]"),
  make_option("--centre_and_scale", action="store_true", type="logical", default=FALSE,
              help="Center and scale the counts by row (genes) [default %default]"),
  make_option("--center_and_scale", action="store_true", type="logical", default=FALSE, dest = "centre_and_scale",
              help="Center and scale the counts by row (genes) [default %default]"),
  make_option("--cluster", action="store", type="character", default="none",
              help="Cluster rows, columns, both or none [default %default]"),
  make_option("--output_clusters", action="store_true", type="logical", default=FALSE,
              help="Output each cluster as a gene list [default %default]"),
  make_option("--cluster_files_base", action="store", type="character", default="clusters-",
              help="Base name for the cluster gene list filenames [default %default]"),
  make_option(c("-k", "--cluster_number"), action="store", type="integer", default=3,
              help="Number of clusters to output (max = 26) [default %default]"),
  make_option("--gene_tree", action="store_true", type="logical", default=FALSE,
              help="Plot dendrogram for gene clustering [default %default]"),
  make_option("--sample_tree", action="store_true", type="logical", default=FALSE,
              help="Plot dendrogram for sample clustering [default %default]"),
  make_option("--colour_palette", type="character", default='plasma',
              help="palette for colour scale [default %default]" ),
  make_option("--cell_colour", type="character", default=NULL,
              help="colour for the outlines of the cells [default %default]" ),
  make_option("--gene_names", action="store_true", type="logical", default=FALSE,
              help="print gene names [default %default]" ),
  make_option("--sample_names", action="store_true", type="logical", default=FALSE,
              help="print sample names [default %default]" ),
  make_option("--sample_metadata_file", type="character", default=NULL,
              help="File of sample metadata to plot below the heatmap [default %default]" ),
  make_option("--sample_metadata_ycol", action="store", default='Category',
              help="Column in sample metadata to plot on y-axis [default %default]"),
  make_option("--sample_metadata_fill", action="store", default='Value',
              help="Column in sample metadata to plot as fill colour [default %default]"),
  make_option("--sample_metadata_fill_palette", action="store", default=NULL,
              help="Fill palette for sample metadata plot [default colour-blind friendly palette from biovisr]"),
  make_option("--gene_metadata_file", type="character", default=NULL,
              help="File of gene metadata to plot to the right of the heatmap [default %default]" ),
  make_option("--gene_metadata_xcol", action="store", default='Category',
              help="Column in gene metadata to plot on x-axis [default %default]"),
  make_option("--gene_metadata_fill", action="store", default='Value',
              help="Column in gene metadata to plot as fill colour [default %default]"),
  make_option("--gene_metadata_fill_palette", action="store", default=NULL,
              help="Fill palette for gene metadata plot [default colour-blind friendly palette from biovisr]"),
  make_option("--relative_widths", action="store", type="character", default="1,3,1",
              help=paste("The relative widths of the plots as a comma-separated list. [default %default]")),
  make_option("--relative_heights", action="store", type="character", default="1,3,1",
              help=paste("The relative heights of the plots as a comma-separated list. [default %default]")),
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
#     "genes_file" = "test_data/test_genes.txt",
#     "detct" = FALSE,
#     "transform" = "rlog",
#     "centre_and_scale" = TRUE,
#     "cluster" = "both",
#     "gene_tree" = TRUE,
#     "sample_tree" = TRUE,
#     "colour_palette" = 'plasma',
#     "cell_colour" = "grey80",
#     "gene_names" = TRUE,
#     "sample_names" = TRUE,
#     "sample_metadata_file" = "test_data/test_samples_long.tsv",
#     "sample_metadata_ycol" = "category",
#     "sample_metadata_fill" = "value",
#     "sample_metadata_fill_palette" = NULL,
#     "gene_metadata_file" = "test_data/test_gene_metadata.tsv",
#     "gene_metadata_xcol" = "category",
#     "gene_metadata_fill" = "value",
#     "gene_metadata_fill_palette" = NULL,
#     "relative_widths" = "1,3,1",
#     "relative_heights" = "1,3,1",
#     "width" = 7,
#     "height" = 10,
#     "fill_limits" = NULL,
#     "output_count_file" = NULL,
#     "output_rda_file" = NULL
#   ),
#   args = c(
#     'test_data/test_samples.tsv',
#     'test_data/test_rnaseq_data.tsv',
#     'test_data/test_heatmap_with_metadata.pdf'
#   )
# )

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'gene_expr_heatmap.R',
    usage = "Usage: %prog [options] sample_file count_file output_file" ),
  positional_arguments = 3
)

output_file <- cmd_line_args$args[3]
packages <- c('tidyverse', 'viridis', 'rnaseqtools', 'biovisr', 'patchwork', 'DESeq2')
if (grepl('svg$', output_file)) { packages <- c(packages, 'svglite') }
metadata_files <- c(cmd_line_args$options[['sample_metadata_file']],
                   cmd_line_args$options[['gene_metadata_file']])
if (any( !sapply(metadata_files, is.null) )) {
    if (any(grepl("\\.ftr$", metadata_files))) { 
      packages <- c(packages, 'feather') 
    }
}
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

fill_palette <- cmd_line_args$options[['colour_palette']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
# unpack relative plot_sizes options
relative_widths <- as.integer(unlist(str_split(cmd_line_args$options[['relative_widths']], ",")))
relative_heights <- as.integer(unlist(str_split(cmd_line_args$options[['relative_heights']], ",")))

if (!is.null(cmd_line_args$options[['fill_limits']])) {
    fill_limits <-
        strsplit(cmd_line_args$options[['fill_limits']], ",", fixed = TRUE)[[1]]
} else {
    fill_limits <- cmd_line_args$options[['fill_limits']]
}

# load sample data
samples <- read_tsv(cmd_line_args$args[1])

#data <- read.delim(cmd_line_args$args[2], check.names = FALSE)
data <- load_rnaseq_data(cmd_line_args$args[2])

# CHECK ALL FILES BEFORE DOING RLOG/VST TRANSFORM
# open genes file if it is specified so that if it doesn't exist the error happens before the rlog/vst transform
if (!is.null(cmd_line_args$options[['genes_file']])) {
    genes <- read_tsv(file = cmd_line_args$options[['genes_file']],
                      col_names = FALSE)
    if (!grepl("ENS[A-Z]*[0-9]+", genes[1,1])){
      # assume there's a header
      genes <- read_tsv(file = cmd_line_args$options[['genes_file']])
    }
}
gene_metadata <- NULL
if (!is.null(cmd_line_args$options[['gene_metadata_file']])) {
  if (grepl("\\.ftr$", cmd_line_args$options[['gene_metadata_file']])) {
    gene_metadata <- read_feather(cmd_line_args$options[['gene_metadata_file']])
  } else {
    gene_metadata <- read_tsv(cmd_line_args$options[['gene_metadata_file']])
  }
}
sample_metadata <- NULL
if (!is.null(cmd_line_args$options[['sample_metadata_file']])) {
  if (grepl("\\.ftr$", cmd_line_args$options[['sample_metadata_file']])) {
    sample_metadata <- read_feather(cmd_line_args$options[['sample_metadata_file']])
  } else {
    sample_metadata <- read_tsv(cmd_line_args$options[['sample_metadata_file']])
  }
}

# make unique rownames depending on DeTCT/RNA-seq 
if (cmd_line_args$options[['detct']]) {
    data$id <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    data$id <- data[['GeneID']]
}

# get count data
counts <- get_counts(data, samples, normalised = FALSE)
# add gene ids to rownames
rownames(counts) <- data$GeneID
# Transform using DESeq2
dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ condition)
dds <- estimateSizeFactors(dds)
if (cmd_line_args$options[['transform']] == "rlog") {
  dds <- rlogTransformation(dds, blind=TRUE)
} else {
  dds <- varianceStabilizingTransformation(dds, blind=TRUE)
}

counts <- assay(dds)
if (cmd_line_args$options[['centre_and_scale']]) {
    counts <- as.data.frame(t( scale( t(counts) ) ))
}

# subset data to gene list if provided
if (!is.null(cmd_line_args$options[['genes_file']])) {
    # subset data to genes and arrange it in same order
    data <- filter(data, GeneID %in% genes[[1]]) %>% 
      arrange(match(GeneID, genes[[1]]))
    # and subset counts
    counts <- counts[ genes[[1]], ]
}

# skip clustering if only 1 row
cluster_rows <- FALSE
cluster_columns <- FALSE
gene_tree <- NULL
sample_tree <- NULL
if (nrow(counts) > 1) {
  # Clustering
  cluster <- function(df) {
      df <- as.data.frame(df)
      distance_matrix <- dist(df)
      clustering <- hclust(distance_matrix)
      df_ordered <- df[ clustering$order, ]
      return(list(df_clustered = df_ordered, clustering_obj = clustering))
  }
  
  if (cmd_line_args$options[['cluster']] == "rows") {
      cluster_rows <- TRUE
  } else if (cmd_line_args$options[['cluster']] == "columns") {
      cluster_columns <- TRUE
  } else if (cmd_line_args$options[['cluster']] == "both") {
      cluster_rows <- TRUE
      cluster_columns <- TRUE
  }
  if (cluster_rows) {
      cluster_list <- cluster(counts)
      counts <- cluster_list[['df_clustered']]
      gene_tree <- cluster_list[['clustering_obj']]
  }
  if (cluster_columns) {
      cluster_list <- cluster(t(counts))
      counts <- as.data.frame(t( cluster_list[['df_clustered']] ))
      sample_tree <- cluster_list[['clustering_obj']]
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
  # make function to return a gene name for an Ensembl id
  ids2names <- function(ids) {
    map_chr(ids, function(id, data){ data$Name[ data$GeneID == id] }, data)
  }
  heatmap_plot <- heatmap_plot + 
    scale_y_discrete(name = NULL, labels = ids2names)
}
heatmap_plot <- heatmap_plot + 
  scale_x_discrete(position = "top") +
  biovisr::theme_heatmap( xaxis_labels = cmd_line_args$options[['sample_names']],
                          yaxis_labels = cmd_line_args$options[['gene_names']],
                          base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust=0)) +
        NULL

# create tree plots if appropriate options are set
gene_tree_plot <- plot_spacer()
if (cluster_rows & cmd_line_args$options[['gene_tree']]) {
  gene_tree_plot <- dendro_plot(gene_tree, categorical_scale = TRUE) + coord_flip() + scale_y_reverse()
}

sample_tree_plot <- plot_spacer()
if (cluster_columns & cmd_line_args$options[['sample_tree']]) {
  sample_tree_plot <- dendro_plot(sample_tree, categorical_scale = TRUE)
}

if (cluster_rows & cmd_line_args$options[['output_clusters']]) {
  groups <- cutree(gene_tree, k = cmd_line_args$options[['cluster_number']])
  clusters <- tibble(GeneID = names(groups), cluster = letters[groups]) %>% 
    mutate(., GeneID = factor(GeneID, levels = levels(plot_data$id))) %>% 
    arrange(., GeneID) %>% 
    mutate(., cluster = factor(cluster, levels = unique(cluster)))
  clusters_list <- split(clusters, clusters$cluster)
  # function to output a gene id file for a cluster
  output_genes <- function(cluster, file_base) {
    filename <- paste0(file_base, cluster$cluster[1], '.txt')
    write.table(cluster$GeneID, file = filename, 
                quote = FALSE, sep = "\t", row.names=FALSE,
                col.names = FALSE)
  }
  lapply(clusters_list, output_genes, cmd_line_args$options[['cluster_files_base']])
  
  if (cmd_line_args$options[['gene_tree']]) {
    gene_tree_plot <- gene_tree_plot + 
      geom_tile(data = clusters, aes(x = GeneID, y = -1, fill = factor(cluster))) +
      scale_fill_manual(values = biovisr::cbf_palette(nlevels(factor(clusters$cluster))),
                        guide = guide_legend(reverse = TRUE),
                        name = "Cluster") + 
      theme(legend.position = "left")
  }
  
}

# load gene metadata if provided
gene_metadata_plot <- plot_spacer()
if (!is.null(gene_metadata)) {
  category_col <- cmd_line_args$options[['gene_metadata_xcol']]
  category_var <- rlang::sym(category_col)
  fill_col <- cmd_line_args$options[['gene_metadata_fill']]
  fill_var <- rlang::sym(fill_col)
  # make xaxis col a factor
  if (class(gene_metadata[[category_col]]) == 'character') {
    gene_metadata[[category_col]] <-
      factor(gene_metadata[[category_col]],
             levels = unique(gene_metadata[[category_col]]))
  }
  
  # set levels of fill column
  gene_metadata[[fill_col]] <-
    factor(gene_metadata[[fill_col]],
           levels = unique(gene_metadata[[fill_col]]))
  
  # subset gene_metadata to genes and order by current gene order
  # and sort by category and then fill
  gene_metadata <- left_join((select(plot_data, id) %>% unique()), 
                             gene_metadata, by = c("id" = "GeneID"))
  # set gene id levels by clustering
  gene_metadata[['id']] <-
    factor(gene_metadata[['id']],
           levels = levels(plot_data$id))
  gene_metadata[[category_col]][ is.na(gene_metadata[[category_col]]) ] <- levels(gene_metadata[[category_col]])[1]
  
  # sort out fill_palette
  fill_palette <- cmd_line_args$options[['gene_metadata_fill_palette']]
  if (!is.null(fill_palette)) {
    if (fill_palette %in% colnames(gene_metadata)) {
      # make a named vector of levels to colours
      if( length(unique(gene_metadata[[fill_palette]])) ==
          nlevels(gene_metadata[[fill_col]]) ) {
        cat2colour <- gene_metadata[ , c(fill_col, fill_palette)] %>% unique() %>% 
          filter(!is.na(!!fill_var))
        fill_palette <- cat2colour[[fill_palette]]
        names(fill_palette) <- cat2colour[[fill_col]]
      }
    }
  }
  
  gene_metadata_plot <-
    df_heatmap(gene_metadata, x = category_col, y = 'id',
               fill = fill_col, fill_palette = fill_palette,
               # colour = "grey50", size = 0.5,
               xaxis_labels = TRUE, yaxis_labels = FALSE,
               na.translate = FALSE
    ) + guides(fill = guide_legend(title = fill_col, reverse = FALSE)) +
    scale_x_discrete(position = "top") +
    theme(axis.title = element_blank(),
          axis.text.x.top = element_text(angle = 45, hjust = 0))
}

# load metadata if provided
# The first column should be the sample ids/names (x axis, matching the samples file)
# the y axis and fill variables are specified by metadata_ycol and metadata_fill
sample_metadata_plot <- plot_spacer()
if (!is.null(sample_metadata)) {
  # subset metadata to samples
  sample_metadata <- filter(sample_metadata, sample %in% levels(plot_data$Sample))
  # order sample levels by clustering
  sample_metadata[['sample']] <-
    factor(sample_metadata[['sample']],
           levels = levels(plot_data$Sample))
  # reverse levels of y axis to make plot match
  category_col <- cmd_line_args$options[['sample_metadata_ycol']]
  if (class(sample_metadata[[category_col]]) == 'character') {
    sample_metadata[[category_col]] <-
      factor(sample_metadata[[category_col]],
              levels = rev(unique(sample_metadata[[category_col]])))
  } else if (class(sample_metadata[[category_col]]) == 'factor') {
    sample_metadata[[category_col]] <-
      factor(sample_metadata[[category_col]],
             levels = rev(levels(sample_metadata[[category_col]])))
  }
  
  # sort out fill_palette
  fill_col <- cmd_line_args$options[['sample_metadata_fill']]
  # set levels of fill column
  sample_metadata[[fill_col]] <-
    factor(sample_metadata[[fill_col]],
            levels = unique(sample_metadata[[fill_col]]))
  
  fill_palette <- cmd_line_args$options[['sample_metadata_fill_palette']]
  if (!is.null(fill_palette)) {
    if (fill_palette %in% colnames(sample_metadata)) {
      # make a named vector of levels to colours
      if( length(unique(sample_metadata[[fill_palette]])) ==
                          nlevels(sample_metadata[[fill_col]]) ) {
        cat2colour <- sample_metadata[ , c(fill_col, fill_palette)] %>% unique()
        fill_palette <- cat2colour[[fill_palette]]
        names(fill_palette) <- cat2colour[[fill_col]]
      }
    }
  }
  
  sample_metadata_plot <-
    df_heatmap(sample_metadata, x = 'sample', y = category_col,
               fill = fill_col, fill_palette = fill_palette,
               colour = "grey50", size = 0.5,
               xaxis_labels = FALSE, yaxis_labels = TRUE,
               na.translate = FALSE
    ) + guides(fill = guide_legend(title = fill_col, reverse = FALSE)) +
    theme(axis.title = element_blank())
}

if (grepl('ps$', output_file)) {
  postscript(file = output_file, paper = "special", horizontal = FALSE,
              width = plot_width, height = plot_height )
} else if (grepl('svg$', output_file)) {
  svglite(file = output_file, width = plot_width, height = plot_height )
} else if (grepl('png$', output_file)) {
  png(file = output_file, width = plot_width, height = plot_height )
} else { # default to pdf 
  pdf(file = output_file, width = plot_width, height = plot_height )
}

plot_list <- list(plot_spacer(), sample_tree_plot, plot_spacer(),
                  gene_tree_plot, heatmap_plot, gene_metadata_plot,
                  plot_spacer(), sample_metadata_plot, plot_spacer())
if (any(class(sample_tree_plot) == "spacer")) {
  relative_heights[1] <- 0
}
if (any(class(gene_tree_plot) == "spacer")) {
  relative_widths[1] <- 0
}
if (any(class(gene_metadata_plot) == "spacer")) {
  relative_widths[3] <- 0
}
if (any(class(sample_metadata_plot) == "spacer")) {
  relative_heights[3] <- 0
}
wrap_plots(
  plot_list,
  ncol = 3,
  nrow = 3,
  byrow = TRUE,
  widths = relative_widths,
  heights = relative_heights,
  guides = 'collect'
)
dev.off()

# reverse counts matrix so it fits the heatmap
counts <- counts[ seq(nrow(counts), 1), ]

if (!is.null(cmd_line_args$options[['output_count_file']])) {
    write.table(counts, file = cmd_line_args$options[['output_count_file']],
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

if (!is.null(cmd_line_args$options[['output_rda_file']])) {
  object_list <- c('data', 'counts', 'plot_data', 'heatmap_plot')
  if (cluster_rows) {
    object_list <- c(object_list, 'gene_tree')
    if (cmd_line_args$options[['gene_tree']]) {
      object_list <- c(object_list, 'gene_tree_plot')
    }
  }
  if (cluster_columns) {
    object_list <- c(object_list, 'sample_tree')
    if (cmd_line_args$options[['sample_tree']]) {
      object_list <- c(object_list, 'sample_tree_plot')
    }
  }
  if (!is.null(cmd_line_args$options[['sample_metadata_file']])) {
    object_list <- c(object_list, 'sample_metadata_plot')
  }
  if (!is.null(cmd_line_args$options[['gene_metadata_file']])) {
    object_list <- c(object_list, 'gene_metadata_plot')
  }
  save(list = object_list, file = cmd_line_args$options[['output_rda_file']])
}
