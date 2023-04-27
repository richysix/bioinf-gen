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
  make_option("--output_file", type="character", default='counts.pdf',
              help="Output file name [default %default]" ),
  make_option("--genes_file", type="character", default=NULL,
              help="File of gene/region ids to subset plots to [default %default]" ),
  make_option("--pvalue_file", type="character", default=NULL,
              help="File of pvalue to add to plots [default %default]" ),
  make_option("--asterisks", action="store_true", type="logical", default=FALSE,
              help="Convert pvalues to asterisks [default %default]" ),
  make_option("--x_variable", type="character", default='condition',
              help="Name of column from samples file to plot on x-axis [default %default]" ),
  make_option("--colour_variable", type="character", default='condition',
              help="Name of column from samples file to plot as colour [default %default]" ),
  make_option("--colour_palette", type="character", default=NULL,
              help=paste("comma-separated list of factor levels to colours (e.g. sib=blue,mut=red or sib=#0000ff,mut=#ff0000) [default %default]",
                          "If NULL the colour-blind friendly palette from biovisr will be used unless there are too many categories",
                          sep = "\n")),
  make_option("--shape_variable", type="character", default=NULL,
              help="Name of column from samples file to plot as shape [default %default]" ),
  make_option("--facet_variable", type="character", default=NULL,
              help="Name of column from samples file to use for facetting [default %default]" ),
  make_option("--crossbar", type="character", default=NULL,
              help="Type of crossbar (mean/median) to add to each group [default %default]" ),
  make_option("--log10", type="logical", action="store_true", default=FALSE,
              help="Use a log10 scaled y-axis [default %default]" ),
  make_option("--width", type="numeric", default=10,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=7,
              help="height of plot (inches) [default %default]" ),
  make_option("--theme_base_size", type="numeric", default=12,
              help="theme_base_size of plot (points) [default %default]" ),
  make_option("--rotate_xaxis_labels", type="logical", action="store_true", default=FALSE,
              help="Rotate x-axis labels to 90 degrees [default %default]" ),
  make_option("--output_data_file", type="character", default=NULL,
              help="Output a Rdata file of the plot object [default %default]" ),
  make_option("--no_jitter", type="logical", action="store_true", default=FALSE,
              help="Don't add jitter to the points  [default %default]" ),
  make_option("--seed", type="integer", default=25673,
              help="random seed to make the jitter reproducible [default %default]" ),
  make_option("--no_pvalue", type="logical", action="store_true", default=FALSE,
              help="Don't add a pvalue to the plots  [default %default]" ),
  make_option("--detct", type="logical", action="store_true", default=FALSE,
              help="Data is DeTCT data, not RNA-seq [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
            help="Turns on debugging statements [default %default]" )
)

# For testing. If running this script interactively the options
# get set to defaults and positional arguments are set to
# whatever is in the arguments vector below
if (any(commandArgs() == "--interactive")) {
  arguments <- 
    file.path('test_data', 
              c('test_samples.tsv', 'test_rnaseq_data.tsv'))
} else {
  arguments <- commandArgs(trailingOnly = TRUE)
}

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'graph_counts_by_group_facet.R',
    usage = "Usage: %prog [options] samples_file count_file" ),
  args = arguments,
  positional_arguments = 2
)
debug <- cmd_line_args$options[['debug']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
theme_base_size <- cmd_line_args$options[['theme_base_size']]
jitter <- !cmd_line_args$options[['no_jitter']]
output_data_file <- cmd_line_args$options[['output_data_file']]
pvalue_file <- cmd_line_args$options[['pvalue_file']]
asterisks <- cmd_line_args$options[['asterisks']]
detct <- cmd_line_args$options[['detct']]

packages <- c('tidyverse', 'biovisr', 'rnaseqtools', 'scales')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# set progress options
options(readr.show_progress = FALSE)

# Read samples
if (debug) { cat("Samples\n") }
samples_file <- cmd_line_args$args[1]
samples <- read.delim( samples_file, header=TRUE, row.names=1 )
# add sample name as factor
samples$sample <- rownames(samples)
samples$sample <- factor(samples$sample,
                        levels = samples$sample)

## NEED TO ADD CHECKING OF COLUMNS. DO SUPPLIED VARIABLES EXIST IN THE DATA

# set levels of colour and shape variable
x_var <- cmd_line_args$options[['x_variable']]
samples[[x_var]] <- factor(samples[[x_var]],
                                levels = unique(samples[[x_var]]))

colour_var <- cmd_line_args$options[['colour_variable']]
samples[[colour_var]] <- factor(samples[[colour_var]],
                                levels = unique(samples[[colour_var]]))

facet_var <- cmd_line_args$options[['facet_variable']]
if (!is.null(facet_var)) {
    samples[[facet_var]] <- factor(samples[[facet_var]],
                                    levels = unique(samples[[facet_var]]))
}
shape_var <- cmd_line_args$options[['shape_variable']]
if (!is.null(shape_var)) {
    samples[[shape_var]] <- factor(samples[[shape_var]],
                                    levels = unique(samples[[shape_var]]))
}

# make colour palette for the data
num_levels <- nlevels(samples[[colour_var]])
if (!is.null(cmd_line_args$options[['colour_palette']])) {
  # split by ',' then by '='
  if (grepl("=", cmd_line_args$options[['colour_palette']])) {
    colours <- str_split_1(cmd_line_args$options[['colour_palette']], ",")  |>  
      str_split(pattern = "=")
    colour_palette <- sapply(colours, function(x){ return(x[2]) })
    names(colour_palette) <- sapply(colours, function(x){ return(x[1]) })
    # check that names match levels
    if (any( !(levels(samples[[colour_var]]) %in% names(colour_palette)) )) {
      missing_levels <- levels(samples[[colour_var]])[ !(levels(samples[[colour_var]]) %in% names(colour_palette)) ]
      print(colour_palette)
      cat(missing_levels, "\n")
      stop('names in colour_palette option do not match levels of colour variable!')
    }
  } else {
    colour_palette <- str_split_1(cmd_line_args$options[['colour_palette']], ",")
  }
  
  # check that number of colours match the number of levels
  if (length(colour_palette) > num_levels) {
    warning('There are more colours than necessary in colour_palette option')
  }
  if (length(colour_palette) < num_levels) {
    stop('Not enough colours in colour_palette option!')
  }
} else if (num_levels <= 10) {
  colour_palette <- cbf_palette(nlevels(samples[[colour_var]]))
  names(colour_palette) <- levels(samples[[colour_var]])
} else{
  ord1 <- seq(1,num_levels,2)
  ord2 <- seq(2,num_levels,2)
  colour_palette <- scales::hue_pal()(num_levels)[ order(c(ord1,ord2)) ]
}
if (debug) {
  print(colour_palette)
}

if (!is.null(shape_var)) {
  make_shape_palette <- function(shape_palette, samples) {
    shape_palette <- shape_palette[seq_len(nlevels(samples[[shape_var]]))]
    names(shape_palette) <- levels(samples[[shape_var]])
    return(shape_palette)
  }
  shape_palette <- 21:25
  if (nlevels(samples[[shape_var]]) <= length(shape_palette)) {
    shape_palette <- make_shape_palette(shape_palette, samples)
  } else {
    shape_palette <- c(15:19,4,6,9)
    if (nlevels(samples[[shape_var]]) <= length(shape_palette)) {
      shape_palette <- make_shape_palette(shape_palette, samples)
    } else {
      stop(paste0('Shape variable (', shape_var, ') has too many levels!'))
    }
  }
}

if (debug) { cat("Counts\n") }

# Read data
data_file <- cmd_line_args$args[2]
data <- load_rnaseq_data(data_file)

# find adjusted pvalue column
if (cmd_line_args$options[['no_pvalue']]) {
  adjp_col <- NULL
} else if (sum(grepl("adjp", names(data))) == 1 ) {
  adjp_col <- names(data)[ grepl("adjp", names(data)) ]
} else {
  adjp_col <- NULL
}

# make a region column
if (detct) {
    data$region <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    data$region <-
        sprintf("%s:%d-%d:%s", data[['Chr']], data[['Start']],
                data[['End']], data[['Strand']])
}

# subset to genes in genes file if it exists
if (!is.null(cmd_line_args$options[['genes_file']])) {
  genes <- read.delim(cmd_line_args$options[['genes_file']], header=FALSE)
  # subset to Gene IDs and arrange in the same order
  data <- filter(data, GeneID %in% genes[[1]]) %>%
    arrange(., match(GeneID, genes[[1]]))
}
regions <- unique(data$region)

# get normalised counts
normalised_counts <- get_counts(data, samples = samples, normalised = TRUE) %>% 
  mutate(GeneID = data$GeneID,
         region = data$region)

# load p value data if it exists
if (!is.null(pvalue_file)) {
  pvalues <- read_tsv(pvalue_file) %>%
    mutate(
      condition1 = factor(condition1, levels = levels(counts_for_plotting$condition)),
      condition2 = factor(condition2, levels = levels(counts_for_plotting$condition)),
      x1 = as.integer(condition1),
      x2 = as.integer(condition2),
      start = case_when(
        x1 < x2 ~ x1,
        x2 < x1 ~ x2
      ),
      midpoint = start + abs(x1 - x2)/2)
  if (asterisks) {
    pvalues <- pvalues %>%
      mutate(
        pval_txt = case_when(
          adjp < 0.001 ~ "***",
          adjp < 0.01 ~ "**",
          adjp < 0.05 ~ "*",
          TRUE ~ "NS")
      )
  } else {
    pvalues <- pvalues %>%
      mutate(
        pval_txt = sprintf("Adjusted p = %.2g, Log2FC = %.2f", adjp, log2fc)
      )
  }
}

make_count_data_for_plot <- function(region_to_plot, normalised_counts, samples){
  plot_data <- normalised_counts |> 
    # subset data to region
    filter(region == region_to_plot) |> 
    # pivot 
    pivot_longer(cols = -c('GeneID', 'region'),
                 names_to = "sample", values_to = "count") |> 
    # set levels of sample
    mutate(sample = factor(sample, levels = unique(sample))) |> 
    # join to samples
    inner_join(samples, by = c('sample'))
  if (cmd_line_args$options[['log10']]) {
    plot_data$log10 <- log10(plot_data$count + 1)
  }
  return(plot_data)
}

# if pdf, open plot device first
output_plot <- function(plot, output_file_name, plot_num) {
  # get output type from filename suffix
  plot_suffix <- sub("^.*\\.", "", output_file_name)
  # pdf is default if nothing matches
  if (plot_suffix == "eps") {
    filename <- sub("\\.eps", paste0('.', plot_num, '.eps'), output_file_name)
    postscript(file = filename, width = plot_width, height = plot_height,
                paper="special", horizontal = FALSE)
    print(plot)
    invisible(dev.off())
  } else if (plot_suffix == "svg") {
    # if svglite is not installed use svg function
    filename <- sub("\\.svg", paste0('.', plot_num, '.svg'), output_file_name)
    if (length(find.package('svglite', quiet = TRUE)) == 0) {
      grDevices::svg(filename = filename, width = plot_width, height = plot_height)
    } else {
      svglite::svglite(file = filename, width = plot_width, height = plot_height)
    }
    print(plot)
    invisible(dev.off())
  } else { # assume pdf device has already been opened and print
    print(plot)
  }
}

make_count_plot <- function(plot_num, data, normalised_counts, samples) {
  region_to_plot <- regions[plot_num]
  if (debug) { cat(region_to_plot, "\n") }
  
  counts <- make_count_data_for_plot(region_to_plot, normalised_counts, samples)
  gene_id <- filter(data, region == region_to_plot) %>% pull("GeneID")
  gene_name <- filter(data, region == region_to_plot) %>% pull("Name")
  if(!is.null(adjp_col)){
    pval <- filter(data, region == region_to_plot) %>% pull(!!adjp_col)
    title <- sprintf("%s\n%s / %s\np = %.3g", region_to_plot,
                     gene_id, gene_name, pval)
  } else {
    title <- sprintf("%s\n%s / %s", region_to_plot,
                     gene_id, gene_name)
  }
  
  # create basic plot
  plot <- ggplot(data = counts, aes_(x = as.name(x_var),
                                     y = as.name('count')))
  
  # add crossbars
  if (!is.null(cmd_line_args$options[['crossbar']])) {
    if (cmd_line_args$options[['crossbar']] == 'mean') {
      plot <- plot + 
        stat_summary(fun = "mean", geom = "crossbar",
                     width = 0.6, size = 0.4 )
    } else if (cmd_line_args$options[['crossbar']] == 'median') {
      plot <- plot + 
        stat_summary(fun = "median", geom = "crossbar",
                     width = 0.6, size = 0.4 )
    }
  }
  
  if (jitter) {
    pos <- position_jitter(width = 0.3, height = 0, seed = cmd_line_args$options[['seed']])
  } else {
    pos <- position_identity()
  }
  # add points with the correct shapes
  if (is.null(shape_var)) {
    plot <- plot +
      geom_point(aes_(fill = as.name(colour_var)),
                 size = 3, shape = 21, colour = 'black',
                 position = pos ) +
      scale_fill_manual(values = colour_palette)
  } else {
    if (shape_palette[1] == 15) {
      plot <- plot +
        geom_point(aes_(colour = as.name(colour_var),
                        shape = as.name(shape_var)), size = 3,
                   position = pos) +
        scale_colour_manual(values = colour_palette,
                            guide = guide_legend(override.aes =
                                                   list(shape = shape_palette[1]),
                                                 order = 1)) +
        scale_shape_manual(values = shape_palette,
                           guide = guide_legend(order = 2))
    } else {
      plot <- plot +
        geom_point(aes_(fill = as.name(colour_var),
                        shape = as.name(shape_var)), size = 3,
                   position = pos ) +
        scale_fill_manual(values = colour_palette,
                          guide = guide_legend(override.aes =
                                                 list(shape = shape_palette[1]),
                                               order = 1)) +
        scale_shape_manual(values = shape_palette,
                           guide = guide_legend(order = 2))
    }
  }
  
  # add pvalue lines and text if necessary
  if (!is.null(pvalue_file)) {
    # subset pvalues to gene/region
    if (detct) {
      pvals <- filter(pvalues, region == region_to_plot)
    } else {
      pvals <- filter(pvalues, GeneID == gene_id)
    }
    # create y values based on largest data value
    points_y_range <- range(counts$count)
    # get major and minor breaks values
    major_breaks <- scales::extended_breaks()(points_y_range)
    minor_breaks <- scales::regular_minor_breaks()(major_breaks, points_y_range, n = 2)
    # set spacing of pvalue lines to halfway between minor breaks
    line_spacing <- (minor_breaks[2] - minor_breaks[1])/2
    n <- nrow(pvals) - 1
    # create y column in pvals by starting at highest major break and adding line spacing
    pvals$ypos <- major_breaks[length(major_breaks)] + seq(0,line_spacing*n,line_spacing)
    
    # with the pval lines added the plot gets bigger so the axes and breaks may change
    # work out minor breaks and how much to nudge the pvalue labels by
    y_limits <- range(c(counts$count, pvals$ypos))
    major_breaks <- scales::extended_breaks()(y_limits)
    minor_breaks <- scales::regular_minor_breaks()(major_breaks, y_limits, n = 2)
    if (asterisks) {
      # nudge asterisk 10% of the distance between the minor breaks
      text_nudge <- (minor_breaks[2] - minor_breaks[1])*0.1
    } else {
      # nudge asterisk 20% of the distance between the minor breaks
      text_nudge <- (minor_breaks[2] - minor_breaks[1])*0.2
    }
    
    # add segments and text to plot
    plot <- plot +
      geom_segment(data = pvals, aes(x = condition1, xend = condition2,
                                     y = ypos, yend = ypos)) +
      geom_text(data = pvals, aes(x = midpoint, y = ypos, label = pval_txt),
                nudge_y = text_nudge)
  }
  
  if (cmd_line_args$options[['log10']]) {
    plot <- plot + scale_y_continuous(trans = scales::pseudo_log_trans()) +
      labs(title = title, y = expression(log[10] * '(Normalised Counts)') )
  } else {
    plot <- plot + labs(title = title, y = "Normalised Counts")
  }
  
  # add facets
  if (!is.null(facet_var)) {
    plot <- plot +
      facet_wrap(as.name(facet_var), nrow = 1)
  }
  
  plot <- plot + 
    theme_minimal(base_size = theme_base_size) +
    theme(strip.background = element_rect(fill = "grey90",
                                          colour = "grey90"),
          axis.text.x = element_text(angle = ifelse(cmd_line_args$options[['rotate_xaxis_labels']], 90, 0)))
  
  output_plot(plot, cmd_line_args$options[['output_file']], plot_num)
}

plot_suffix <- sub("^.*\\.", "", cmd_line_args$options[['output_file']])
if(plot_suffix == 'pdf') {
  pdf(file = cmd_line_args$options[['output_file']],
      width = plot_width, height = plot_height)
}

purrr::walk(seq_along(regions), make_count_plot, data, normalised_counts, samples)

if(plot_suffix == 'pdf') {
  dev.off()
}

# # save rds file if output_data_file options is set
# if (!is.null(output_data_file)) {
#   save(plot_list, samples, data, counts_for_plotting, file = output_data_file)
# }
