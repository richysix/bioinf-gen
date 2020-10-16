library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='go_bubble.pdf',
              help="Output file name [default %default]" ),
  make_option("--label_topn", type="integer", default=5,
              help="Label the top n points (by pvalue) [default %default]" ),
  make_option("--label_p_cutoff", type="numeric", default=NULL,
              help="Label points with pvalue below the cutoff [default %default]" ),
  make_option("--labels", type="character", default=NULL,
              help="comma-separated list of GO.IDs indicating which points to label [default %default]" ),
  make_option("--no_labels", action="store_true", type="logical", default=FALSE,
              help="comma-separated list of GO.IDs indicating which points to label [default %default]" ),
  make_option("--width", type="numeric", default=7,
              help="width of plot (inches) [default %default]" ),
  make_option("--height", type="numeric", default=10,
              help="height of plot (inches) [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to turn topgo output into a bubble chart',
  paste('The script assumes there are files named',
        '"BP.sig.tsv", "CC.sig.tsv" and "MF.sig.tsv"',
        'in the current working directory'),
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'go_bubble_plot.R',
    usage = "Usage: %prog [options] ",
    description = desc ),
  positional_arguments = 0
)

# cmd_line_args <-
#   list(options = list('output_file' = 'go_bubble.pdf',
#                       'label_topn' = 5,
#                       'labels' = 'GO:0016208,GO:0006955,GO:0098844',
#                       'label_p_cutoff' = 0.0001,
#                       'no_labels' = FALSE,
#                       'width' = 10,
#                       'height' = 7,
#                       debug = TRUE),
#        args = c())

# load packages
packages <- c('tidyverse', 'viridis', 'biovisr', 'miscr', 'ggrepel')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

go_bp <- read_tsv('BP.sig.tsv') %>% 
  mutate(., domain = "BP")
go_mf <- read_tsv('MF.sig.tsv') %>% 
  mutate(., domain = "MF")
go_cc <- read_tsv('CC.sig.tsv') %>% 
  mutate(., domain = "CC")
go_results <- rbind(go_mf, go_bp, go_cc) %>% 
  mutate(., log10p = -log10(pval),
         label = domain)

if (cmd_line_args$options[['no_labels']]) {
  labels_l <- FALSE
} else if (!is.null(cmd_line_args$options[['labels']])) {
  labels_l <- TRUE
  go_terms_to_label <- unlist(str_split(cmd_line_args$options[['labels']], ","))
  go_results <- go_results %>% 
    mutate(label = case_when(GO.ID %in% go_terms_to_label ~ "label",
                             TRUE ~ as.character(label)))
} else if (!is.null(cmd_line_args$options[['label_p_cutoff']])) {
  labels_l <- TRUE
  go_results <- go_results %>% 
    mutate(label = case_when(pval <= cmd_line_args$options[['label_p_cutoff']] ~ "label",
                             TRUE ~ as.character(label)))
} else {
  labels_l <- TRUE
  go_terms_to_label <- go_results %>% 
    arrange(., desc(log10p)) %>% 
    head(., cmd_line_args$options[['label_topn']]) %>% 
    pull(., var = GO.ID)
  go_results <- go_results %>% 
    mutate(label = case_when(GO.ID %in% go_terms_to_label ~ "label",
                             TRUE ~ as.character(label)))
}
if (cmd_line_args$options[['debug']]) {
  print(go_results %>% filter(., label == "label") %>% pull(., GO.ID))
}

# go_bubble_plot <- ggplot(data = go_results, aes(x = GO.ID, y = log10p)) +
#   geom_point(aes(size = Significant, colour = domain),
#              alpha = 0.8) +
#   geom_text_repel(data = go_results[ go_results$label == "label", ],
#                   aes(label = Term), size = 8/72*25.4,
#                   nudge_x = 0.2, nudge_y = 0.2,
#                   hjust = 0, direction = "y") +
#   facet_wrap(vars(domain), nrow = 1,
#              strip.position = "bottom") +
#   scale_x_discrete(name = NULL, expand = expand_scale(mult = 0.05)) +
#   scale_color_manual(values = biovisr::cbf_palette(length(unique(go_results$domain))),
#                      guide  = "none") +
#   scale_size_area(guide = guide_legend(override.aes = list(alpha = 0.6)),
#                         name = "Number of\nSignificant Genes") + 
#   labs(y = expression(-log[10]*"[pvalue]")) +
#   theme_minimal() + 
#   theme(axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank())

# TEST highlight labelled points as well
fill_palette <- biovisr::cbf_palette(length(unique(go_results$domain)))
names(fill_palette) <- sort(unique(go_results$domain))
if (labels_l) {
  colour_palette <- c(fill_palette, "label" = "black")
} else {
  colour_palette <- fill_palette
}

go_bubble_plot_highlighted <- ggplot(data = go_results, aes(x = GO.ID, y = log10p)) +
  geom_point(aes(size = Significant, fill = domain, colour = label),
             alpha = 0.8, shape = 21) +
  geom_text_repel(data = go_results[ go_results$label == "label", ],
                  aes(label = Term), size = 10/72*25.4,
                  nudge_x = 0.2, nudge_y = 0.2,
                  hjust = 0, direction = "y") +
  facet_wrap(vars(domain), nrow = 1,
             strip.position = "bottom") +
  scale_x_discrete(name = NULL, expand = expand_scale(mult = 0.05)) +
  scale_color_manual(values = colour_palette,
                     guide  = "none") +
  scale_fill_manual(values = fill_palette,
                     guide  = "none") +
  scale_size_area(guide = guide_legend(override.aes = list(alpha = 0.6, shape = 16)),
                  name = "Number of\nSignificant Genes") +
  labs(y = expression(-log[10]*"[pvalue]")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank())
# png(filename = 'go_bubble_highlight_test.png',
#     width = 960, height = 540)
# print(go_bubble_plot_highlighted)
# dev.off()

# print plot to file
output_plot(list(plot = go_bubble_plot, 
                 filename = cmd_line_args$options[['output_file']]),
            width = cmd_line_args$options[['width']],
            height = cmd_line_args$options[['height']])
