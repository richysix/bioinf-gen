#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='bubble_plot.pdf',
              help="Name of output file [default %default]" ),
  make_option("--x_var", type="character", default='x',
              help="Column to plot on the x axis [default %default]" ),
  make_option("--x_levels", type="character", default=NULL,
              help="Ordering for x axis (comma-separated list) [default the order they appear in the input file]" ),
  make_option("--y_var", type="character", default='y',
              help="Column to plot on the y axis [default %default]" ),
  make_option("--size_var", type="character", default='size',
              help="Column to map to size [default %default]" ),
  make_option("--fill_var", type="character", default='fill',
              help="Column to map to fill [default %default]" ),
  make_option("--x_labels", type="character", default=NULL,
              help="Column to use as labels for the x axis [default %default]" ),
  make_option("--y_labels", type="character", default=NULL,
              help="Column to use as labels for the y axis [default %default]" ),
  make_option("--reverse_y", action="store_true", default=FALSE,
              help="Switch to reverse the ordering of the y axis [default %default]" ),
  make_option("--width", type="integer", default=NULL,
              help="Width of the plot [default %default]" ),
  make_option("--height", type="integer", default=NULL,
              help="Height of the plot [default %default]" ),
  make_option("--output_data_file", type="character", default=NULL,
              help="Output a Rdata file of the plot object [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'bubble_plot.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

# convert options to less cumbersome variable names
x_col <- cmd_line_args$options[['x_var']]
y_col <- cmd_line_args$options[['y_var']]
size_col <- cmd_line_args$options[['size_var']]
fill_col <- cmd_line_args$options[['fill_var']]
x_lab_col <- cmd_line_args$options[['x_labels']]
y_lab_col <- cmd_line_args$options[['y_labels']]
reverse_y <- cmd_line_args$options[['reverse_y']]
output_file <- cmd_line_args$options[['output_file']]
width <- cmd_line_args$options[['width']]
height <- cmd_line_args$options[['height']]
output_data_file <- cmd_line_args$options[['output_data_file']]

if (!is.null(cmd_line_args$options[['x_levels']])) {
  x_levels <- unlist(strsplit(cmd_line_args$options[['x_levels']],
                              ",", fixed = TRUE))
} else {
  x_levels <- NULL
}

packages <- c('tidyverse', 'ggplot2', 'biovisr')
# add svglite to packages if output file name ends in svg
if (grepl("svg$", output_file)) {
  packages <- c(packages, 'svglite')
}
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# load input data
input_data <- read_tsv(cmd_line_args$args[1])
# set levels of input_data
# check class of x and y data
# if character, convert to factors in the order that they appear in the file
# set levels of input_data
# check class of x and y data
# if character, convert to factors in the order that they appear in the file
# unless x_levels is set
if( class(input_data[[ x_col ]]) == 'character' ) {
  if (!is.null(x_levels)) {
    # TO DO: CHECK X_LEVELS MATCHES LEVELS IN DATA
    input_data[[ x_col ]] <- factor(input_data[[ x_col ]],
                                    levels = x_levels)
  } else {
    input_data[[ x_col ]] <- factor(input_data[[ x_col ]],
                                    levels = unique(input_data[[ x_col ]]))
  }
}

# reverse the order for the y axis if option --reverse_y is TRUE
if( class(input_data[[ y_col ]]) == 'character' ) {
  if( reverse_y ) {
    y_levels <- rev(unique(input_data[[ y_col ]]))
  } else {
    y_levels <- unique(input_data[[ y_col ]])    
  }
  input_data[[ y_col ]] <- factor(input_data[[ y_col ]],
                                levels = y_levels)
}

# check labels
if ( is.null(x_lab_col) ) {
  x_labels <- NULL 
} else {
  x_levels <- levels(input_data[[ x_col ]])
  x_lab_df <- unique(input_data[, c(x_col, x_lab_col) ])
  x_labels <- x_lab_df[[x_lab_col]]
  names(x_labels) <- x_lab_df[[x_col]]
}

# Need to make sure the y's and y labels are the same length
if ( is.null(y_lab_col) ) {
  y_labels <- NULL 
} else {
  y_lab_df <- unique(input_data[, c(y_col, y_lab_col) ])
  if( reverse_y ) {
    y_labels <- rev(y_lab_df[[ y_lab_col ]])
  } else {
    y_labels <- y_lab_df[[ y_lab_col ]]
  }
}

plot <- 
  biovisr::bubble_plot(input_data, x = x_col, y = y_col,
                    size = size_col, fill = fill_col,
                    x_labels = x_labels,
                    y_labels = y_labels)

# plot to file
if (grepl("png$", output_file)) {
  width <- ifelse(is.null(width), 480, width)
  height <- ifelse(is.null(height), 480, height)
  png(file = output_file, width = width, height = height)
} else if (grepl("ps$", output_file)) {
  width <- ifelse(is.null(width), 10, width)
  height <- ifelse(is.null(height), 11, height)
  postscript(file = output_file, width = width, height = height,
             paper = "special", horizontal = FALSE)
} else if (grepl("svg$", output_file)) {
  width <- ifelse(is.null(width), 10, width)
  height <- ifelse(is.null(height), 11, height)
  svglite(file = output_file, width = width, height = height)
} else {
  width <- ifelse(is.null(width), 10, width)
  height <- ifelse(is.null(height), 11, height)
  pdf(file = output_file, width = width, height = height)
}
print(plot)
dev.off()

# save rds file if output_data_file options is set
if (!is.null(output_data_file)) {
  save(plot, input_data, file = output_data_file)
}
