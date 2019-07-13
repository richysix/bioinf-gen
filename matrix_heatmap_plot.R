#!/usr/bin/env Rscript

library(optparse)

# parse options
option_list <- list(
  make_option(c("-o", "--output_file"), type="character", default='heatmap.pdf',
              help="Output File name [default %default]" ),
  make_option("--x_column", type="character", default='1',
              help="column number/name to use as the x axis variable [default %default]" ),
  make_option("--x_labels_column", type="character", default=NULL,
              help="column to use as axis tick labels for x axis [default %default]" ),
  make_option("--y_column", type="character", default='2',
              help="column number/name to use as the y axis variable [default %default]" ),
  make_option("--y_labels_column", type="character", default=NULL,
              help="column to use as axis tick labels for y axis [default %default]" ),
  make_option("--data_column", type="character", default='6',
              help="column number/name to plot as colour [default %default]" ),
  make_option("--colour_palette", type="character", default=NULL,
              help="palette for colour scale [default %default]" ),
  make_option("--colour_legend_position", type="character", default='right',
              help="position for colour scale legend [default %default]" ),
  make_option("--na_colour", type="character", default='grey95',
              help="colour for NA values [default %default]" ),
  make_option("--reverse_x_axis", action="store_true", default=FALSE,
              help="Reverse the order of the x axis [default %default]" ),
  make_option("--reverse_y_axis", action="store_true", default=FALSE,
              help="Reverse the order of the y axis [default %default]" ),
  make_option("--x_axis_label", type="character", default=NULL,
              help="label for x axis [default %default]" ),
  make_option("--x_axis_position", type="character", default='bottom',
              help="position of x axis ('bottom' or 'top')[default %default]" ),
  make_option("--y_axis_label", type="character", default=NULL,
              help="label for y axis [default %default]" ),
  make_option("--y_axis_position", type="character", default='left',
              help="position of x axis ('left' or 'right')[default %default]" ),
  make_option("--data_axis_label", type="character", default=NULL,
              help="label for colour scale [default %default]" ),
  make_option("--header", action="store_true", default=FALSE,
              help="Does the input file have a header line [default %default]" ),
  make_option("--width", type="numeric", default=7,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=10,
              help="height of plot (inches) [default %default]" ),
  make_option("--output_data_file", type="character", default=NULL,
              help="Output a Rdata file of the plot object [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'matrix_heatmap_plot.R',
    usage = "Usage: %prog [options] inputFile",
    description = "Generic heatmap script. The input file is expected to be in long form. \nIt requires at a minimum a column for the x variable, one for the y variable and a data column" ),
  positional_arguments = 1
)

# load packages
packages <- c('ggplot2')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# only load viridis if scale is selected
if( !is.null(cmd_line_args$options[['colour_palette']]) ) {
  if (cmd_line_args$options[['colour_palette']] == 'viridis'){
    suppressPackageStartupMessages(
      suppressWarnings( library('viridis', character.only = TRUE) )
    )
  }
}

#print(cmd_line_args)

if( cmd_line_args$options[['verbose']] ){
  cat( "Output File:", cmd_line_args$options[['outputFile']], "\n", sep=" " )
  cat( "Number/Name of X column:", cmd_line_args$options[['x_column']], "\n", sep=" " )
  cat( "Number/Name of Y column:", cmd_line_args$options[['y_column']], "\n", sep=" " )
  cat( "Number/Name of data column:", cmd_line_args$options[['data_column']], "\n", sep=" " )
  cat( "Colour Palette:", cmd_line_args$options[['colour_palette']], "\n", sep=" " )
  cat( "Reverse X axis:", cmd_line_args$options[['reverse_x_axis']], "\n", sep=" " )
  cat( "Reverse Y axis:", cmd_line_args$options[['reverse_y_axis']], "\n", sep=" " )
  cat( "Input has header:", cmd_line_args$options[['header']], "\n", sep=" " )
  cat( "Label for X axis:", cmd_line_args$options[['x_axis_label']], "\n", sep=" " )
  cat( "Label for Y axis:", cmd_line_args$options[['y_axis_label']], "\n", sep=" " )
  cat( "Label for data scale:", cmd_line_args$options[['data_axis_label']], "\n", sep=" " )
}

# load data
input_data <- read.table(file = cmd_line_args$args[1], sep = "\t",
                         header = cmd_line_args$options[['header']])
# if header exists match column names by name, else match by number
if ( cmd_line_args$options[['header']] ) {
  x_column <- cmd_line_args$options[['x_column']]
  y_column <- cmd_line_args$options[['y_column']]
  data_column <- cmd_line_args$options[['data_column']]
  if (!is.null(cmd_line_args$options[['x_labels_column']])) {
    x_labels_column <- cmd_line_args$options[['x_labels_column']]
  } else {
    x_labels_column <- NULL
  }

  if (!is.null(cmd_line_args$options[['y_labels_column']])) {
    y_labels_column <- cmd_line_args$options[['y_labels_column']]
  } else {
    y_labels_column <- NULL
  }
} else {
  x_column <- names(input_data)[ as.integer(cmd_line_args$options[['x_column']]) ]
  y_column <- names(input_data)[ as.integer(cmd_line_args$options[['y_column']]) ]
  data_column <- names(input_data)[ as.integer(cmd_line_args$options[['data_column']]) ]
  if (!is.null(cmd_line_args$options[['x_labels_column']])) {
    x_labels_column <- names(input_data)[ as.integer(cmd_line_args$options[['x_labels_column']]) ]
  } else {
    x_labels_column <- NULL
  }
  
  if (!is.null(cmd_line_args$options[['y_labels_column']])) {
    y_labels_column <- names(input_data)[ as.integer(cmd_line_args$options[['y_labels_column']]) ]
  } else {
    y_labels_column <- NULL
  }
}

# set levels of x and y as order in which they appear in input file
# then reverse if reverse_x_axis option set
input_data[[x_column]] <-
  factor(input_data[[x_column]], levels = unique(input_data[[x_column]]))
if ( cmd_line_args$options[['reverse_x_axis']] ) {
  input_data[[x_column]] <-
    factor(input_data[[x_column]],
      levels = rev(levels(input_data[[x_column]])))
}
if (!is.null(x_labels_column)) {
  input_data[[x_labels_column]] <-
    factor(input_data[[x_labels_column]],
      levels = unique(input_data[[x_labels_column]]))
  if ( cmd_line_args$options[['reverse_x_axis']] ) {
    input_data[[x_labels_column]] <-
      factor(input_data[[x_labels_column]],
        levels = rev(levels(input_data[[x_labels_column]])))
  }
  x_labels <- levels(input_data[[x_labels_column]])
} else {
  x_labels <- NULL
}

# reverse levels of x-axis if option set
input_data[[y_column]] <-
  factor(input_data[[y_column]], levels = unique(input_data[[y_column]]))
if ( cmd_line_args$options[['reverse_y_axis']] ) {
  input_data[[y_column]] <-
    factor(input_data[[y_column]],
      levels = rev(levels(input_data[[y_column]])))
}
if (!is.null(y_labels_column)) {
  input_data[[y_labels_column]] <-
    factor(input_data[[y_labels_column]],
      levels = unique(input_data[[y_labels_column]]))
  if ( cmd_line_args$options[['reverse_y_axis']] ) {
    input_data[[y_labels_column]] <-
      factor(input_data[[y_labels_column]],
        levels = rev(levels(input_data[[y_labels_column]])))
  }
  y_labels <- levels(input_data[[y_labels_column]])
} else {
  y_labels <- NULL
}

# make plot
heatmap_plot <- ggplot(data = input_data) +
  geom_raster( aes_q(x = as.name(x_column), y = as.name(y_column),
                     fill = as.name(data_column) ) ) +
  theme_void() +
  theme(axis.title = element_text(size = 14, colour = 'black', angle = 0, debug = FALSE),
        axis.text.x = element_text(size = 12, colour = 'black', angle = 45,
                                   hjust = 0, vjust = 0.5, debug = FALSE),
        axis.text.y = element_text(size = 12, colour = 'black', angle = 0, debug = FALSE),
        legend.title = element_text(size = 14, colour = 'black', angle = 0, debug = FALSE),
        legend.text = element_text(size = 12, colour = 'black', angle = 0, debug = FALSE),
        legend.position = cmd_line_args$options[['colour_legend_position']])

if( !is.null(cmd_line_args$options[['colour_palette']]) ) {
  if(cmd_line_args$options[['colour_palette']] == 'viridis'){
    heatmap_plot <- heatmap_plot + scale_fill_viridis(na.value = cmd_line_args$options[['na_colour']])
  } else if (cmd_line_args$options[['colour_palette']] == 'diverging') {
    heatmap_plot <- heatmap_plot +
      scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b',
                           midpoint = 0,
                           na.value = cmd_line_args$options[['na_colour']])
  }
}
if (!is.null(x_labels) ) {
  heatmap_plot <- heatmap_plot +
    scale_x_discrete(position =  cmd_line_args$options[['x_axis_position']],
                     labels = x_labels)
} else {
  heatmap_plot <- heatmap_plot +
    scale_x_discrete(position =  cmd_line_args$options[['x_axis_position']])
}
if (!is.null(y_labels) ) {
  heatmap_plot <- heatmap_plot +
    scale_y_discrete(position =  cmd_line_args$options[['y_axis_position']],
                     labels = y_labels)
} else {
  heatmap_plot <- heatmap_plot +
    scale_y_discrete(position =  cmd_line_args$options[['y_axis_position']])
}

if ( !is.null(cmd_line_args$options[['x_axis_label']]) ) {
  heatmap_plot <- heatmap_plot + xlab(cmd_line_args$options[['x_axis_label']])
}
if ( !is.null(cmd_line_args$options[['y_axis_label']]) ) {
  heatmap_plot <- heatmap_plot + ylab(cmd_line_args$options[['y_axis_label']])
}
if ( !is.null(cmd_line_args$options[['data_axis_label']]) ) {
  heatmap_plot <- heatmap_plot + guides( fill = guide_colourbar(title = cmd_line_args$options[['data_axis_label']] ) )
}

if(grepl('pdf$', cmd_line_args$options[['output_file']])) {
  pdf(file = cmd_line_args$options[['output_file']],
      width = cmd_line_args$options[['width']],
      height = cmd_line_args$options[['height']] )
} else if (grepl('ps$', cmd_line_args$options[['output_file']])) {
  postscript(file = cmd_line_args$options[['output_file']],
              paper = "A4", horizontal = FALSE,
              width = cmd_line_args$options[['width']],
              height = cmd_line_args$options[['height']] )
} else if (grepl('svg$', cmd_line_args$options[['output_file']])) {
  library(svglite)
  svglite(file = cmd_line_args$options[['output_file']],
          width = cmd_line_args$options[['width']],
          height = cmd_line_args$options[['height']] )
} else { # default to pdf 
  pdf(file = cmd_line_args$options[['output_file']],
      width = cmd_line_args$options[['width']],
      height = cmd_line_args$options[['height']] )
}
print(heatmap_plot)
dev.off()

if (!is.null(cmd_line_args$options[['output_data_file']])) {
  save(heatmap_plot, file = cmd_line_args$options[['output_data_file']])
}
