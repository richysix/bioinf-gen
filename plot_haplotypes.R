#!/usr/bin/env Rscript

packages <- c("ggplot2", "scales", "optparse", "svglite")
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

option_list <- list(
  make_option("--basename", type="character", default='haplotypes',
              help="Base of output files [default %default]"),
  make_option("--haplotype_levels", type="character", default='DHTu2,HET,DHAB',
              help="Comma separated list of the levels of the haplotype column in order [default %default]"),
  make_option("--haplotype_colours", type="character", default='blue,green,red',
              help="Comma separated list of the colours for the levels of haplotype [default %default]"),
  make_option("--orientation", type="character", default='landscape',
              help="Orientation of chromosome (landscape = chr on x-axis) [default %default]"),
  make_option("--output_type", type="character", default='pdf',
              help="Type of file to output [default %default]"),
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output")
)

# ADD PLOT WIDTH AND HEIGHT AS OPTIONS

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'plot_haplotypes.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

# cmd_line_args <- list(
#   options = list(directory="/Users/rjw26/citiid/hpc-work/sat-ase/call_haplotypes",
#                  haplotype_levels = "DHTu2/DHTu2,DHAB/DHTu2,DHAB/DHAB,NC,NA",
#                  haplotype_colours = "#0073B3,#009980,#CC6600,grey50,grey80",
#                  orientation = "portrait",
#                  basename='sat_rnaseq-zmp_ph209-sat_5_6-haplotypes-1Mb', verbose=TRUE ),
#   args = c("~/citiid/hpc-work/sat-ase/call_haplotypes/sat_rnaseq-zmp_ph209-sat_5_6-haplotypes-1Mb-long.txt") )

if (!is.null(cmd_line_args$options[['directory']])) {
  output_base <- file.path(cmd_line_args$options[['directory']], cmd_line_args$options[['basename']])
} else {
  output_base <- file.path(cmd_line_args$options[['basename']])
}

# OPTIONS
if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
  cat( "Output files basename:", cmd_line_args$options[['basename']], "\n", sep=" " )
  cat( "Haplotype levels:", cmd_line_args$options[['haplotype_levels']], "\n", sep=" " )
  cat( "Haplotype colours:", cmd_line_args$options[['haplotype_colours']], "\n", sep=" " )
}

# haplotype levels
haplotype_levels <- unlist(strsplit(cmd_line_args$options[['haplotype_levels']], ","))
# haplotype_colours
haplotype_colours <- unlist(strsplit(cmd_line_args$options[['haplotype_colours']], ","))
names(haplotype_colours) <- haplotype_levels

if ("NA" %in% names(haplotype_colours)) {
  na_colour <- haplotype_colours["NA"]
  haplotype_colours <- haplotype_colours[ setdiff(names(haplotype_colours), c("NA")) ]
} else {
  na_colour <- "grey80"
}

# orientation
orientation <- cmd_line_args$options[['orientation']]
if (!(orientation %in% c('landscape', 'portrait'))) {
  stop("--orientation must be one of landscape or portrait")
}

# read in data
input_file <- cmd_line_args$args[1]
if (input_file == "-") {
  input_data <- read.delim(file = file("stdin"))
} else {
  input_data <- read.delim(file = input_file)
}

# order haplotype
input_data$haplotype <- factor(input_data$haplotype, levels = haplotype_levels)
# order sample
input_data$sample <- factor(input_data$sample,
                            levels = rev(unique(input_data$sample)))

# set up centres
input_data$x <- input_data$start + (input_data$end - input_data$start)/2

# widths
input_data$w <- input_data$end - input_data$start

to_Mb <- function(vec){
  vec <- as.character( as.numeric(vec) / 1000000 )
}
plot_haplotypes <- function( plot_data, orientation = "landscape",
                             haplotype_colours = biovisr::cbf_palette(),
                             na_colour = 'grey80'){
  x_lab = paste0("Position (Mb) - Chr ", plot_data$chr[1])
  if (orientation == "portrait") {
    plot_data$sample <- factor(plot_data$sample, levels = rev(levels(plot_data$sample)))
  }
  hap_plot <- ggplot(data = plot_data) + 
    geom_tile( aes( x = x, y = sample, width = w, fill = haplotype ) ) + 
    scale_fill_manual(values = haplotype_colours, na.value = na_colour,
                      name = "Haplotype") + 
    labs(x=x_lab, y="Sample") +
    scale_x_continuous(labels = to_Mb) +
    theme_void()
  if (orientation == "landscape") {
    hap_plot <- hap_plot + 
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.title.y = element_text(angle = 90))
  } else {
    hap_plot <- hap_plot + 
      coord_flip() +
      theme(axis.text = element_text(size = 16),
            # axis.text.x = element_text(size = 16, angle = 45, 
            #                            vjust = 1, hjust = 1),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 18),
            axis.title.y = element_text(angle = 90),
            # legend.position = "top",
            # legend.direction = "vertical",
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16))
  }
  return(hap_plot)
}

plot_list <- lapply( split(input_data, input_data$chr), plot_haplotypes,
                     orientation = cmd_line_args$options[["orientation"]],
                     haplotype_colours = haplotype_colours,
                     na_colour = na_colour)

# output plots to file
if (cmd_line_args$options[['output_type']] == 'svg') {
  for (i in seq_len(length(plot_list))) {
    svg_filename <- paste(paste0(output_base, '-Chr', i), 'svg', sep = '.')
    svglite(file = svg_filename, width=7, height=5)
    print(plot_list[[i]])
    dev.off()
  }
} else if (cmd_line_args$options[['output_type']] == 'eps') {
  for (i in seq_len(length(plot_list))) {
    eps_filename <- paste(paste0(output_base, '-Chr', i), 'eps', sep = '.')
    postscript(file = eps_filename, width=7, height=5,
               horizontal = TRUE, paper = 'special')
    print(plot_list[[i]])
    dev.off()
  }
} else { # use pdf as fallback
  pdf_filename <- paste(output_base, "pdf", sep=".")
  pdf(file=pdf_filename, width=9, height=10, paper = "special")
  invisible(print(plot_list))
  invisible(dev.off())
}
