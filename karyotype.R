library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default = 'ideogram.png',
              help="Output file name [default %default]" ),
  make_option("--heatmap_data", type="character", default=NULL,
              help="tsv of data for chromosome heatmap. Expects columns Chr, Start, End, Value [default %default]" ),
  make_option("--heatmap_colours", type="character", default=NULL,
              help="comma separated string of colours for heatmap [default %default]" ),
  make_option("--annotation_data", type="character", default='mutation.txt',
              help="tsv of annotation data [default %default]" ),
  make_option("--annotation_label_type", type="character", default='marker',
              help="Type of annotation track (marker, line, polygon, heatmap) [default %default]" ),
  make_option("--width", type="integer", default=NULL,
              help="Width of the plot [default %default]" ),
  make_option("--height", type="integer", default=NULL,
              help="Height of the plot [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to plot a karyotype ideogram',
  sep = "\n"
)

if (interactive()) {
  cmd_args <- c("/Users/rjw26/work/apocrita/data/scratch/bty114/zmp-network/Dr.e85.chrom.sizes")
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}

cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list, prog = 'karyotype.R',
    description = desc,
    usage = "Usage: %prog [options] chr.sizes"),
  args = cmd_args,
  positional_arguments = 1
)

# load packages
packages <- c('tidyverse', 'RIdeogram', 'biovisr')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load chromosome sizes file
chr_sizes <- read_tsv(cmd_line_args$args[1],
                      col_names = c("Chr", "Start", "End"))

# load heatmap data
if (is.null(cmd_line_args$options[['heatmap_data']])) {
  heatmap_data <- NULL
} else {
  heatmap_data <- read_tsv(cmd_line_args$options[['heatmap_data']])
}
# load colours
if (is.null(cmd_line_args$options[['heatmap_colours']])) {
  heatmap_colours <- c("yellow", "orange", "red")
} else {
  heatmap_colours <- str_split(cmd_line_args$options[['heatmap_colours']], ",")[[1]]
}

# load annotation data
if (is.null(cmd_line_args$options[['annotation_data']])) {
  annotation_data <- NULL
} else {
  annotation_data <- read_tsv(cmd_line_args$options[['annotation_data']])
}

# plot ideogram
plot_device <- sub("^.*\\.", "", cmd_line_args$options[['output_file']])
output_svg <- sub("\\..*$", ".svg", cmd_line_args$options[['output_file']])
plot_width <- ifelse(is.null(cmd_line_args$options[['width']]), 8.2677, cmd_line_args$options[['width']])
plot_height <- ifelse(is.null(cmd_line_args$options[['height']]), 11.6929, cmd_line_args$options[['height']])

ideogram(karyotype = chr_sizes, output = output_svg,
         # overlaid = heatmap_data, colorset1 = heatmap_colours,
         label = annotation_data, label_type = cmd_line_args$options[['annotation_label_type']])
convertSVG(output_svg, device = plot_device,
           file = cmd_line_args$options[['output_file']],
           width = plot_width, height = plot_height )
