#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='go_barchart.pdf',
              help="Output file name [default %default]" ),
  make_option("--x_variable", type="character", default='FE',
              help="Name of column to plot on the x axis [default %default]" ),
  make_option("--x_axis_title", type="character", default='Fold Enrichment',
              help="Title for the x axis [default %default]" ),
  make_option("--fill_variable", type="character", default='Category',
              help="Name of column to plot as the bar fill colour [default %default]" ),
  make_option("--up_down_variable", type="character", default='up_down',
              help="Name of column to use to determine direction of bars [default %default]" ),
  make_option("--up_down_levels", type="character", default='Down=left,Up=right',
              help="Levels of the up down column with directions [default %default]" ),
  make_option("--width", type="numeric", default=NULL,
              help="Width of the plot [default: 480px for .png, 8 inches for other formats]" ),
  make_option("--height", type="numeric", default=NULL,
              help="Height of the plot [default: 360px for .png, 6 inches for other formats]" ),
  make_option("--top_terms", type="integer", default=NULL,
              help="How many terms to plot (ordered by x variable) [default All terms]" ),
  make_option("--RData_file", type="character", default=NULL,
              help="Name of RData file containing the plots to output [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to plot the GO enrichment results as a bar chart',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'go_barchart.R',
    usage = "Usage: %prog [options] input_files",
    description = desc ),
  positional_arguments = 1
)

# cmd_line_args <-
#   list(options = list('output_file' = 'go_plot.pdf',
#                       'x_variable' = 'FE',
#                       'fill_variable' = 'Set',
#                       'up_down_variable' = 'up_down',
#                       'up_down_levels' = 'Down=left,Up=right',
#                       'top_terms' = NULL,
#                       'RData_file' = 'go_plot.rda',
#                       debug = TRUE),
#        args = c('test_data_go.tsv'))

# load packages
packages <- c('tidyverse', 'rlang', 'biovisr', 'miscr', 'grid')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}
# if(!suppressWarnings(require('extrafont'))) {
#   extrafont <- FALSE
#   font_name <- "Helvetica"
# } else {
#   extrafont <- TRUE
#   font_name <- "Arial"
# }

# set up plot dimensions
output_file <- cmd_line_args$options[['output_file']]
plot_width <- cmd_line_args$options[['width']]
plot_height <- cmd_line_args$options[['height']]
file_suffix <- sub("^.*\\.", "", output_file)
if (file_suffix == "png") {
  plot_width <- ifelse(is.null(plot_width), 480, plot_width)
  plot_height <- ifelse(is.null(plot_height), 360, plot_height)
} else {
  plot_width <- ifelse(is.null(plot_width), 8, plot_width)
  plot_height <- ifelse(is.null(plot_height), 6, plot_height)
}

# load data
go_info <- read_tsv(cmd_line_args$args[1])

# If top_terms is NULL make it equal to the number of terms
if (is.null(cmd_line_args$options[['top_terms']])) {
  num_top_terms <- nrow(go_info)
} else {
  num_top_terms <- cmd_line_args$options[['top_terms']]
}

# go through options
up_down_levels_l <- strsplit(unlist(strsplit(cmd_line_args$options[['up_down_levels']], ",")), "=")
x_var <- rlang::sym(cmd_line_args$options[['x_variable']])
up_down_var <- rlang::sym(cmd_line_args$options[['up_down_variable']])
fill_var <- rlang::sym(cmd_line_args$options[['fill_variable']])

for (up_down_level in up_down_levels_l) {
  if (up_down_level[2] == 'left') {
    go_info <- go_info %>% 
      mutate(., !!x_var := case_when(!!up_down_var == up_down_level[1] ~ -!!x_var,
                               TRUE ~ !!x_var))
  }
}

# arrange by absolute FE/pval
# take the top n terms
# rearrange by FE/pval
# and make GO.ID, up_down_var and fill_var into factors
go_info <- go_info %>% 
  arrange(., desc(abs(!!x_var))) %>% 
  head(., num_top_terms) %>% 
  arrange(., !!x_var) %>% 
  mutate(GO.ID = factor(GO.ID, levels = unique(GO.ID)),
         up_down = factor(!!up_down_var),
         !!fill_var := factor(!!fill_var,
                              levels = unique(go_info[[ cmd_line_args$options[['fill_variable']] ]]) ) )

# Text for labelling the x axis
# up_label <- grobTree(textGrob(label = "Up genes", gp = gpar(fontsize = 12, fontfamily = font_name)))
up_label <- grobTree(textGrob(label = "Up genes", gp = gpar(fontsize = 12)))
# down_label <- grobTree(textGrob(label = "Down genes", gp = gpar(fontsize = 12, fontfamily = font_name)))
down_label <- grobTree(textGrob(label = "Down genes", gp = gpar(fontsize = 12)))

# vectors for labelling the bars
up_term_labels <- go_info %>% 
  filter(., !!x_var > 0) %>% 
  select(., GO.ID, Term) %>% 
  unique(.) %>% 
  mutate(., y = -1)
down_term_labels <- go_info %>% 
  filter(., !!x_var < 0) %>% 
  select(., GO.ID, Term) %>% 
  unique(.) %>% 
  mutate(., y = 1)

# make bar chart
coloured_bar_chart <- ggplot(data = go_info) +
  geom_col(aes(x = GO.ID, y = !!x_var, fill = !!fill_var),
           position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = up_term_labels, hjust = 1, vjust = 0.5, size = 12/72*25.4,
            aes(x = GO.ID, y = y, label = Term)) +
  geom_text(data = down_term_labels, hjust = 0, vjust = 0.5, size = 12/72*25.4,
            aes(x = GO.ID, y = y, label = Term)) +
  scale_fill_manual(values = cbf_palette(nlevels(go_info[[ cmd_line_args$options[['fill_variable']] ]]))) +
  scale_y_continuous(labels = abs) +
  labs(y = cmd_line_args$options[['x_axis_title']], x = "GO Term") +
  annotation_custom(up_label,
                    xmin = -1, xmax = -3,
                    ymin = 10, ymax = 20 ) +
  annotation_custom(down_label,
                    xmin = -1, xmax = -3,
                    ymin = -15, ymax = -5 ) +
  coord_flip(clip = "off") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "top", 
        # text = element_text(family = font_name),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        panel.grid.major.y = element_blank())

# output plot
miscr::open_graphics_device(output_file,
                            width = plot_width,
                            height = plot_height)
print(coloured_bar_chart)
dev.off()

# embed fonts if it's a pdf
# if (grepl("pdf$", cmd_line_args$options[['output_file']])) {
#   if(extrafont) {
#     embed_fonts(cmd_line_args$options[['output_file']], 
#                 outfile = cmd_line_args$options[['output_file']])
#   }
# }

# save to rdata file if option specified
if (!is.null(cmd_line_args$options[['RData_file']])) {
  save(coloured_bar_chart, go_info,
       file = cmd_line_args$options[['RData_file']])
}
