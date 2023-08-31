#!/usr/bin/env Rscript

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

library('optparse')

option_list <- list(
  make_option("--output", type="character", default="all-cor.tsv",
              help="Name of output file [default %default]"),
  make_option("--type", type="character", default="genes",
              help="Make network from genes or samples [default %default]"),
  make_option("--correlation_measure", type="character", default="spearman",
              help="Correlation coefficent to use [default %default]"),
  make_option("--threshold", type="numeric", default=NULL,
              help="Threshold for correlation coefficent to include in network [default %default]"),
  make_option("--clusters_file", type="character", default="all-cor-clusters.tsv",
              help="Name of cluster output file from MCL clustering [default %default]"),
  make_option("--expansion", type="numeric", default=2,
              help="Value of MCL expansion parameter [default %default]"),
  make_option("--inflation", type="numeric", default=1.4,
              help="Value of MCL inflation parameter [default %default]")
)

desc <- paste('', 'Script to create a correlation network from RNA-seq count data', 
              'and cluster it with MCL',
              sep = "\n")
if (interactive()) {
  cmd_args <- c("samples.tsv", "test-first-100.tsv")
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'create_cor_network.R',
    description = desc,
    usage = "Usage: %prog [options] samples_file counts_file" ),
  args = cmd_args,
  positional_arguments = 2
)

cmd_line_args$options[['correlation_measure']] <- tolower(cmd_line_args$options[['correlation_measure']])
if (!(cmd_line_args$options[['correlation_measure']] %in% c('pearson', 'spearman'))) {
  stop('Option --correlation_measure needs to be one of "pearson" and "spearman"')
}
cmd_line_args$options[['type']] <- tolower(cmd_line_args$options[['type']])
if (!(cmd_line_args$options[['type']] %in% c('genes', 'samples'))) {
  stop('Option --type needs to be one of "genes" and "samples"')
}

packages <- c('tidyverse', 'rnaseqtools', 'matrixStats', 'MCL')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load RNAseq data
sample_info <- read_tsv(cmd_line_args$args[1],
                        col_names = c("sample", "condition"),
                        show_col_types = FALSE)

rnaseq_data <- load_rnaseq_data(cmd_line_args$args[2])

# remove any genes where counts are all zero
std_dev <- get_counts(rnaseq_data, samples = sample_info, normalised = TRUE) |> 
  as.matrix() |> 
  rowSds()
rnaseq_data <- rnaseq_data |> 
  filter(std_dev > 0)

# get counts
norm_counts <- get_counts(rnaseq_data, samples = sample_info, normalised = TRUE)

# calculate cor coefficient
if (cmd_line_args$options[['type']] == "samples") {
  cor_m <- cor(norm_counts, method = cmd_line_args$options[['correlation_measure']])
  rownames(cor_m) <- sample_info$sample
  colnames(cor_m) <- sample_info$sample
} else {
  cor_m <- norm_counts |> 
    t() |> 
    cor(method = cmd_line_args$options[['correlation_measure']]) |> 
    abs()
  if (!is.null(cmd_line_args$options[['threshold']])) {
    for (i in seq_len(nrow(cor_m))) {
      for (j in seq_len(ncol(cor_m))) {
        if (cor_m[i,j] < cmd_line_args$options[['threshold']]) {
          cor_m[i,j] <- 0
        }
      }
    }
  }
  rownames(cor_m) <- rnaseq_data$GeneID
  colnames(cor_m) <- rnaseq_data$GeneID
}

# output correlations as edges file
cor_m |> 
  as_tibble(rownames = "gene1") |> 
  pivot_longer(cols = -gene1, names_to = "gene2", values_to = "cor") |> 
  filter(cor > cmd_line_args$options[['threshold']]) |> 
  write_tsv(file = cmd_line_args$options[['output']])

# cluster network with MCL
cor_m_clustered <- mcl(cor_m, expansion = cmd_line_args$options[['expansion']],
    inflation = cmd_line_args$options[['inflation']],
    addLoops = TRUE)

# output clusters file and edges file
clusters <- tibble(
  Cluster_id = as.character(cor_m_clustered$Cluster),
  GeneID = rownames(cor_m)
)
select(rnaseq_data, GeneID:Description) |> 
  inner_join(clusters, by = "GeneID") |> 
  mutate(Cluster_id = fct_infreq(Cluster_id)) |> 
  arrange(Cluster_id) |> 
  write_tsv(file = cmd_line_args$options[['clusters_file']])
