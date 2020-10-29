library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='go_plot.pdf',
              help="Output file name [default %default]" ),
  make_option("--RData_file", type="character", default='go_plot.rda',
              help="Name of RData file containing the plots to output [default %default]" ),
  make_option("--top_terms", type="integer", default=NULL,
              help="How many terms to plot (ordered by Fold Enrichment) [default All terms]" ),
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
  positional_arguments = c(1, Inf)
)

# cmd_line_args <-
#   list(options = list('output_file' = 'go_plot.pdf',
#                       'RData_file' = 'go_plot.rda',
#                       'top_terms' = NULL,
#                       debug = TRUE),
#        args = c('BP.sig.tsv', 'CC.sig.tsv', 'MF.sig.tsv'))

# load packages
packages <- c('tidyverse', 'biovisr', 'grid', 'extrafont')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load data
go_info <- read_tsv(cmd_line_args$args[1])
# make column, combo of Response_Category and Term_ID
go_info <- go_info %>% 
  mutate(., UID = paste(Response_Category, Term_ID, sep ="-"),
         Response_Category = factor(Response_Category, levels = unique(Response_Category))) %>% 
  arrange(., UID)

# If top_terms is NULL make it equal to the number of terms
if (is.null(cmd_line_args$options[['top_terms']])) {
  num_top_terms <- nrow(go_info)
} else {
  num_top_terms <- cmd_line_args$options[['top_terms']]
}

run_fishers_test <- function(line) {
  in_term_in_set <- line$Num_genes_in_term_in_set[1]
  in_term_not_set <- line$Num_genes_in_term[1] - in_term_in_set
  not_term_in_set <- line$Study_set_size[1] - in_term_in_set
  not_term_not_set <- line$Reference_set_size[1] - in_term_in_set - in_term_not_set - not_term_in_set
  fisher.test(cbind(c(in_term_in_set, in_term_not_set), c(not_term_in_set, not_term_not_set)),
              alternative = "greater")
}

fishers_results <- split(go_info, go_info$UID) %>% 
  lapply(., run_fishers_test) %>% 
  sapply(., function(x){ signif(x$p.value, digits = 6) })

if (any(signif(go_info$pvalue, digits = 6) != fishers_results)) {
  ids <- go_info$UID[ (signif(go_info$pvalue, digits = 6) != fishers_results) ]
  go_info[ (signif(go_info$pvalue, digits = 6) != fishers_results), ] %>% 
    cbind(., fisher_results = fishers_results[ids]) %>% 
    select(., UID, pvalue, fisher_results, adjp, Term_name) %>% 
    mutate(., diff = pvalue - fisher_results) %>% 
    print()
  warning("Some of the Pvalues don't match!")
}

go_info$FE <- (go_info$Num_genes_in_term_in_set / go_info$Study_set_size) /
  (go_info$Num_genes_in_term/go_info$Reference_set_size)
# order by fold enrichment and filter to top terms
go_info_top <- go_info %>% 
  arrange(., FE) %>% 
  head(., num_top_terms) %>% 
  mutate(Term_name = factor(Term_name, levels = unique(Term_name)),
         cluster = factor(cluster))

# output plots to pdf
pdf(file = cmd_line_args$options[['output_file']],
    width = 10)
# facetted barchart
facetted_barchart <- ggplot(data = go_info_top) +
  geom_col(aes(x = Term_name, y = FE, fill = cluster)) +
  facet_wrap(vars(Response_Category)) +
  scale_fill_manual(name = NULL,
                    values = cbf_palette(nlevels(go_info_top$cluster)),
                    labels = c('Up', 'Down')) +
  coord_flip() + 
  labs(y = "Fold Enrichment", x = "GO Term") +
  theme(legend.position = "top",
        text = element_text(family = "Arial"))
print(facetted_barchart)

# Coloured bar chart
# change FE to -ve if cluster2
go_info_top <- go_info_top %>% 
  mutate(., FE = ifelse(cluster == 'Cluster1', FE, -FE)) %>% 
  arrange(., FE) %>% 
  mutate(Term_name = factor(Term_name, levels = unique(Term_name)))

up_label <- grobTree(textGrob(label = "Up genes", gp = gpar(fontsize = 10, fontfamily = "Arial")))
down_label <- grobTree(textGrob(label = "Down genes", gp = gpar(fontsize = 10, fontfamily = "Arial")))

# version with labels next to bars
up_term_labels <- go_info_top %>% 
  filter(., FE > 0) %>% 
  select(., Term_name) %>% 
  unique(.) %>% 
  mutate(., y = -1)
down_term_labels <- go_info_top %>% 
  filter(., FE < 0) %>% 
  select(., Term_name) %>% 
  unique(.) %>% 
  mutate(., y = 1)
coloured_bar_chart <- ggplot(data = go_info_top) +
  geom_col(aes(x = Term_name, y = FE, fill = Response_Category),
           position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = up_term_labels, hjust = 1, vjust = 0.5, size = 10/72*25.4,
            aes(x = Term_name, y = y, label = Term_name)) +
  geom_text(data = down_term_labels, hjust = 0, vjust = 0.5, size = 10/72*25.4,
            aes(x = Term_name, y = y, label = Term_name)) +
  scale_fill_manual(values = cbf_palette(nlevels(go_info_top$Response_Category)),
                    ) +
  coord_flip(clip = "off") +
  annotation_custom(up_label,
                    xmin = 0, xmax = -2,
                    ymin = 10, ymax = 20 ) +
  annotation_custom(down_label,
                    xmin = 0, xmax = -2,
                    ymin = -15, ymax = -5 ) +
  labs(y = "Fold Enrichment", x = "GO Term") +
  theme_minimal() +
  theme(legend.position = "top", 
        text = element_text(family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        panel.grid.major.y = element_blank())
print(coloured_bar_chart)

dev.off()

# embed fonts
embed_fonts(cmd_line_args$options[['output_file']], 
            outfile = cmd_line_args$options[['output_file']])

save(facetted_barchart, coloured_bar_chart, go_info_top,
     file = cmd_line_args$options[['RData_file']])
