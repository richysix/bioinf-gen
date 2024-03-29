#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='counts.pdf',
              help="Output file name [default %default]" ),
  make_option("--x_variable", type="character", default='condition',
              help="Name of column from samples file to plot on x-axis [default %default]" ),
  make_option("--group_variable", type="character", default='condition',
              help="Name of column from samples file to group observations by [default %default]" ),
  make_option("--colour_palette", type="character", default=NULL,
              help="comma-separated list of factor levels to colours (e.g. sib=blue,mut=red or sib=#0000ff,mut=#ff0000) [default %default]" ),
  make_option("--errorbars", type="logical", action="store_true", default=FALSE,
              help="Add error bars [default %default]" ),
  make_option("--width", type="numeric", default=10,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=7,
              help="height of plot (inches) [default %default]" ),
  make_option("--theme_base_size", type="numeric", default=12,
              help="theme_base_size of plot (inches) [default %default]" ),
  make_option("--detct", type="logical", action="store_true", default=FALSE,
              help="Data is DeTCT data, not RNA-seq [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
            help="Turns on debugging statements [default %default]" )
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'graph_counts_line.R',
    usage = "Usage: %prog [options] samples_file count_file" ),
  positional_arguments = 2
)
debug <- cmd_line_args$options[['debug']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
theme_base_size <- cmd_line_args$options[['theme_base_size']]

packages <- c('tidyverse', 'biovisr')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# Read samples
if (debug) { cat("Samples\n") }
samples_file <- cmd_line_args$args[1]
samples <- read.table( samples_file, header=TRUE, row.names=1 )
# set levels of colour and shape variable
x_var <- cmd_line_args$options[['x_variable']]
samples[[x_var]] <- factor(samples[[x_var]],
                                levels = unique(samples[[x_var]]))

group_var <- cmd_line_args$options[['group_variable']]
samples[[group_var]] <- factor(samples[[group_var]],
                                levels = unique(samples[[group_var]]))

# add sample name as factor
samples$sample <- rownames(samples)
samples$sample <- factor(samples$sample,
                        levels = samples$sample)

if (debug) { cat("Counts\n") }

# Read data
data_file <- cmd_line_args$args[2]
data <- load_rnaseq_data(data_file)
# data <- read.delim(data_file, header=TRUE, check.names=FALSE)

# # Support different column names
# names(data)[names(data) == 'chr']               <- 'Chr'
# names(data)[names(data) == '#Chr']              <- 'Chr'
# names(data)[names(data) == 'start']             <- 'Start'
# names(data)[names(data) == 'end']               <- 'End'
# names(data)[names(data) == 'strand']            <- 'Strand'
# names(data)[names(data) == 'ID']                <- 'Gene ID'
# names(data)[names(data) == 'adjpval']           <- 'adjp'
# names(data)[names(data) == 'padj']              <- 'adjp'
# names(data)[names(data) == 'Adjusted p value']  <- 'adjp'
# names(data)[names(data) == 'Gene name']         <- 'Name'
# names(data)[ grepl("e[0-9]+ Ensembl Gene ID", names(data)) ] <- 'Gene ID'
# names(data)[names(data) == 'GeneID']                <- 'Gene ID'

# find adjusted pvalue column
if (sum(grepl("adjp", names(data))) == 1) {
  adjp_col <- names(data)[ grepl("adjp", names(data)) ] 
} else {
  adjp_col <- NULL
}

# make a region column
if (cmd_line_args$options[['detct']]) {
    data$region <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    data$region <-
        sprintf("%s:%d-%d:%s", data[['Chr']], data[['Start']],
                data[['End']], data[['Strand']])
}
regions <- unique(data$region)

normalised_counts <- data[ , grepl("region|normalised", colnames(data)) ]
colnames(normalised_counts) <- sub(".normalised.count", "", colnames(normalised_counts))

# order counts by sample
normalised_counts <- normalised_counts[ , c('region', rownames(samples)) ]

# melt counts and join with sample info
if (debug) { cat("Join\n") }
counts_for_plotting <-
    gather(normalised_counts, key = 'sample', value = 'count', -region)
counts_for_plotting$sample <- factor(counts_for_plotting$sample,
                                    levels = unique(counts_for_plotting$sample))
counts_for_plotting <- merge(samples, counts_for_plotting)

# summarise data for lines
calculate_error_for_ci <- function(values, conf = 0.95) {
  p <- 1 - (1 - conf)/2
  qnorm(p)*sd(values)/sqrt(length(values))
}

count_summary <- 
  group_by(counts_for_plotting, stage, condition, region) %>%
  summarise(mean_count = mean(count),
            error = calculate_error_for_ci(count)
          ) %>%
  mutate(lower_ci = mean_count - error,
         upper_ci = mean_count + error)

# set up colour palette
colour_blind_palette <- 
  c( 'blue' = rgb(0,0.45,0.7),
     'vermillion' = rgb(0.8, 0.4, 0),
     'blue_green' = rgb(0, 0.6, 0.5),
     'yellow' = rgb(0.95, 0.9, 0.25),
     'sky_blue' = rgb(0.35, 0.7, 0.9),
     'orange' = rgb(0.9, 0.6, 0),
     'purple' = rgb(0.8, 0.6, 0.7),
     'black' = rgb(0, 0, 0) )
# make colour palette for the data
num_levels <- nlevels(samples[[group_var]])
if (!is.null(cmd_line_args$options[['colour_palette']])) {
    # split by ',' then by '='
    if (grepl("=", cmd_line_args$options[['colour_palette']])) {
        colours <- strsplit(unlist(strsplit(cmd_line_args$options[['colour_palette']], ",")), "=")
        colour_palette <- sapply(colours, function(x){ return(x[2]) })
        names(colour_palette) <- sapply(colours, function(x){ return(x[1]) })
        
        # check that names match levels
        missing_levels <- unlist( lapply(levels(samples[[group_var]]),
                                        function(level, colour_palette){
                                            if( !(level %in% names(colour_palette)) ){
                                                return(level)
                                            }
                                        }, colour_palette ) )
        if (!is.null(missing_levels)) {
            print(colour_palette)
            cat(missing_levels, "\n")
            stop('names in colour_palette option do not match levels of colour variable!')
        }
    } else {
        colour_palette <- unlist(strsplit(cmd_line_args$options[['colour_palette']], ","))
    }
    
    # check that number of colours match the number of levels
    if (length(colour_palette) < num_levels) {
        stop('Not enough colours in colour_palette option!')
    }
} else if (nlevels(samples[[group_var]]) <= length(colour_blind_palette)) {
    colour_palette <- colour_blind_palette[seq_len(nlevels(samples[[group_var]]))]
    names(colour_palette) <- levels(samples[[group_var]])
} else{
    num_levels <- nlevels(samples[[group_var]])
    ord1 <- seq(1,num_levels,2)
    ord2 <- seq(2,num_levels,2)
    colour_palette <- hue_pal()(num_levels)[ order(c(ord1,ord2)) ]
}
if (debug) {
    print(colour_palette)
}

# function to produce a line graph
line_plot <- function(counts, x_var, group_var, colour_palette, errorbars = FALSE) {
  line_size <- 1.5
  plot <- 
    ggplot(data = counts, 
         aes_(x = as.name(x_var), y = as.name('mean_count'), 
             colour = as.name(group_var))) + 
    geom_line(aes_(group = as.name(group_var)), size = line_size) + 
    scale_colour_manual(values = colour_palette)
  if (errorbars) {
    plot <- plot + 
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                    width = 0.2, size = line_size)
  }
  return(plot)
}

# make a plot for each region
plot_list <- lapply(regions,
    function(region, counts_for_plotting, data) {
        if (debug) { cat(region, "\n") }
        counts <- count_summary[ count_summary$region == region, ]
        if (is.null(adjp_col)) {
          title <- sprintf("%s\n%s / %s", region,
                           data[ data$region == region, "Gene ID"],
                           data[ data$region == region, "Name"])
        } else {
          title <- sprintf("%s\n%s / %s\np = %.3g", region,
                           data[ data$region == region, "Gene ID"],
                           data[ data$region == region, "Name"],
                           data[ data$region == region, adjp_col ])
        }
        
        plot <- line_plot(counts, x_var, group_var, 
                          colour_palette, errorbars = cmd_line_args$options[['errorbars']])
        
        plot <- plot +      
            labs(title = title, y = "Normalised Counts") +
            theme_minimal(base_size = theme_base_size)
        return(plot)
    }, counts_for_plotting, data
)

# get output type from filename suffix
# pdf is default if nothing matches
if (sub("^.*\\.", "", cmd_line_args$options[['output_file']]) == "eps") {
    invisible(
        lapply(seq_len(length(plot_list)),
            function(i, plot_list) {
                filename <- sub("\\.eps", paste0('.', i, '.eps'),
                                cmd_line_args$options[['output_file']])
                postscript(file = filename, width = plot_width, height = plot_height,
                           paper="special", horizontal = FALSE)
                print(plot_list[[i]])
                dev.off()
                return(TRUE)
            }, plot_list)
    )
} else if (sub("^.*\\.", "", cmd_line_args$options[['output_file']]) == "svg") {
    library('svglite')
    invisible(
        lapply(seq_len(length(plot_list)),
            function(i, plot_list) {
                filename <- sub("\\.svg", paste0('.', i, '.svg'),
                                cmd_line_args$options[['output_file']])
                svglite(file = filename, width = plot_width, height = plot_height)
                print(plot_list[[i]])
                dev.off()
                return(TRUE)
            }, plot_list)
    )    
} else {
    pdf(file = cmd_line_args$options[['output_file']],
        width = plot_width, height = plot_height)
    invisible(lapply(plot_list, print))
    dev.off()
}
