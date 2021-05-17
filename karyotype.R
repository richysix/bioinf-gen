library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='ideogram.png',
              help="Output file name [default %default]" ),
  make_option("--heatmap_data", type="character", default=NULL,
              help="tsv of data for chromosome heatmap. Expects columns Chr, Start, End, Value [default %default]" ),
  make_option("--heatmap_colours", type="character", default=NULL,
              help="comma separated string of colours [default %default]" ),
  make_option("--annotation_file", type="character", default='mutation.txt',
              help="tsv of annotation [default %default]" ),
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

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'karyotype.R',
    usage = "Usage: %prog [options] chr.sizes",
    description = desc ),
  positional_arguments = 1
)

# cmd_line_args <-
#   list(options = list('output_file' = 'ideogram.png',
#                       'heatmap_data' = 'binned_sig_genes.bed',
#                       debug = TRUE),
#        args = c('~/citiid/hpc-work/resources/genomes/danio_rerio/grcz11/GRCz11.bed'))
 
# load packages
packages <- c('tidyverse', 'RIdeogram')
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
  heatmap_colours <- c("#FFFFFF", "firebrick3")
} else {
  heatmap_colours <- str_split(cmd_line_args$options[['heatmap_colours']], ",")[[1]]
}


# # load annotation data
# annotation <- read_tsv(cmd_line_args$options[['annotation_file']])

# plot ideogram
plot_device <- sub("^.*\\.", "", cmd_line_args$options[['output_file']])
output_svg <- sub("\\..*$", ".svg", cmd_line_args$options[['output_file']])
plot_width <- ifelse(is.null(cmd_line_args$options[['width']]), 8.2677, cmd_line_args$options[['width']])
plot_height <- ifelse(is.null(cmd_line_args$options[['height']]), 11.6929, cmd_line_args$options[['height']])

ideogram(karyotype = chr_sizes, output = output_svg,
         overlaid = heatmap_data, colorset1 = heatmap_colours)
convertSVG(output_svg, device = plot_device,
           file = cmd_line_args$options[['output_file']],
           width = plot_width, height = plot_height )

# , 
#          overlaid = binned_gene_counts, colorset1 = ,
#          label = annotation, label_type = "marker")

