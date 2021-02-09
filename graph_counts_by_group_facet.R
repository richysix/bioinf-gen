#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='counts.pdf',
              help="Output file name [default %default]" ),
  make_option("--genes_file", type="character", default=NULL,
              help="File of gene/region ids to subset plots to [default %default]" ),
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
              help="theme_base_size of plot (inches) [default %default]" ),
  make_option("--rotate_xaxis_labels", type="logical", action="store_true", default=FALSE,
              help="Rotate x-axis labels to 90 degrees [default %default]" ),
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

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'graph_counts_by_group_facet.R',
    usage = "Usage: %prog [options] samples_file count_file" ),
  positional_arguments = 2
)
debug <- cmd_line_args$options[['debug']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
theme_base_size <- cmd_line_args$options[['theme_base_size']]
jitter <- !cmd_line_args$options[['no_jitter']]

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
if (cmd_line_args$options[['detct']]) {
    data$region <- paste(data[['Chr']], data[['Region start']],
                         data[['Region end']], data[["3' end position"]],
                         data[["3' end strand"]], sep=":")
} else {
    data$region <-
        sprintf("%s:%d-%d:%s", data[['Chr']], data[['Start']],
                data[['End']], data[['Strand']])
}

if (!is.null(cmd_line_args$options[['genes_file']])) {
  genes <- read.delim(cmd_line_args$options[['genes_file']], header=FALSE)
  # subset to Gene IDs and arrange in the same order
  data <- filter(data, GeneID %in% genes[[1]]) %>%
    arrange(., match(GeneID, genes[[1]]))
}
regions <- unique(data$region)

# get normalised counts
normalised_counts <- get_counts(data, samples = samples, normalised = TRUE) %>% 
  mutate(region = data$region)

if (debug) { cat("Join\n") }
# subset data to regions if necessary
# normalised_counts <- filter(normalised_counts, region %in% regions)
# melt counts and join with sample info
counts_for_plotting <-
    gather(normalised_counts, key = 'sample', value = 'count', -region)
counts_for_plotting$sample <- factor(counts_for_plotting$sample,
                                    levels = unique(counts_for_plotting$sample))
counts_for_plotting <- merge(samples, counts_for_plotting)

if (cmd_line_args$options[['log10']]) {
  counts_for_plotting$log10 <- log10(counts_for_plotting$count + 1)
}

# make colour palette for the data
num_levels <- nlevels(samples[[colour_var]])
if (!is.null(cmd_line_args$options[['colour_palette']])) {
    # split by ',' then by '='
    if (grepl("=", cmd_line_args$options[['colour_palette']])) {
        colours <- strsplit(unlist(strsplit(cmd_line_args$options[['colour_palette']], ",")), "=")
        colour_palette <- sapply(colours, function(x){ return(x[2]) })
        names(colour_palette) <- sapply(colours, function(x){ return(x[1]) })
        
        # check that names match levels
        missing_levels <- unlist( lapply(levels(samples[[colour_var]]),
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
    if (length(colour_palette) > num_levels) {
      warning('There are more colours than necessary in colour_palette option')
    }
    if (length(colour_palette) < num_levels) {
        stop('Not enough colours in colour_palette option!')
    }
} else if (nlevels(samples[[colour_var]]) <= 10) {
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

plot_list <- lapply(regions,
    function(region, counts_for_plotting, data) {
        if (debug) { cat(region, "\n") }
        counts <- counts_for_plotting[ counts_for_plotting$region == region, ]
        if(!is.null(adjp_col)){
          title <- sprintf("%s\n%s / %s\np = %.3g", region,
                           data[ data$region == region, "GeneID"],
                           data[ data$region == region, "Name"],
                           data[ data$region == region, adjp_col ])
        } else {
          title <- sprintf("%s\n%s / %s", region,
                           data[ data$region == region, "GeneID"],
                           data[ data$region == region, "Name"])
        }
        
        plot <- ggplot(data = counts, aes_(x = as.name(x_var),
                                            y = as.name('count')))

        # add crossbars
        if (!is.null(cmd_line_args$options[['crossbar']])) {
            if (cmd_line_args$options[['crossbar']] == 'mean') {
                plot <- plot + 
                    stat_summary(fun.y = "mean", fun.ymin = "mean",
                                 fun.ymax="mean", geom = "crossbar",
                                 width = 0.6, size = 0.4 )
            } else if (cmd_line_args$options[['crossbar']] == 'median') {
                plot <- plot + 
                    stat_summary(fun.y = "median", fun.ymin = "median",
                                 fun.ymax="median", geom = "crossbar",
                                 width = 0.6, size = 0.4 )
            }
        }
        
        if (jitter) {
          pos <- position_jitter(width = 0.3, height = 0, seed = cmd_line_args$options[['seed']])
        } else {
          pos <- position_identity()
        }
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
