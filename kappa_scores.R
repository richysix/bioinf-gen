#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='kappa_scores.tsv',
              help="Output file name [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to calculate kappa scores from topgo output',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'kappa_scores.R',
    usage = "Usage: %prog [options] sig_file topgo_sig_genes_file",
    description = desc ),
  positional_arguments = 2
)

# cmd_line_args <-
#   list(options = list('output_file' = 'kappa_scores.tsv',
#                       debug = TRUE),
#        args = c('sig.tsv', 'BP.sig.genes.tsv'))

# load packages
packages <- c('tidyverse')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load sig genes
sig_genes <- read_tsv(cmd_line_args$args[1])

# topgo sig.genes.tsv
go_sig_genes <- read_tsv(cmd_line_args$args[2])

go_terms_matrix <- 
  left_join(sig_genes, go_sig_genes, by = c('Gene')) %>% 
  mutate(., GO.ID = case_when(is.na(GO.ID) ~ go_sig_genes$GO.ID[1], 
                                TRUE ~ GO.ID),
           pvalue = case_when(is.na(`p value`) ~ FALSE, 
                              `p value` == 1 ~ TRUE)) %>% 
  select(., Gene, GO.ID, pvalue) %>% 
  pivot_wider(., names_from = GO.ID, values_from = pvalue,
              values_fill = list(pvalue = FALSE))

calc_kappa_score <- function(first_col, second_col, go_terms_matrix){
  num_genes <- nrow(go_terms_matrix)
  InBoth <- sum(go_terms_matrix[[first_col]] & go_terms_matrix[[second_col]])
  InFirstNotSecond <- sum(go_terms_matrix[[first_col]] & !go_terms_matrix[[second_col]])
  InSecondNotFirst <- sum(!go_terms_matrix[[first_col]] & go_terms_matrix[[second_col]])
  InNeither <- sum(!go_terms_matrix[[first_col]] & !go_terms_matrix[[second_col]])
  po <- (InBoth + InNeither)/ num_genes
  go1 <- ((InBoth + InFirstNotSecond)/num_genes)*((InBoth + InSecondNotFirst)/num_genes)
  go2 <- ((InSecondNotFirst + InNeither)/num_genes)*((InFirstNotSecond + InNeither)/num_genes)
  pe <- go1 + go2
  kappa <- (po - pe) / (1 - pe)
  return(kappa)
}

term_ids <- select(go_sig_genes, GO.ID) %>% arrange() %>% unique() %>% pull()
num_terms <- length(term_ids)
pairwise_terms_1 <- character(length = num_terms * (num_terms - 1)/2)
pairwise_terms_2 <- character(length = num_terms * (num_terms - 1)/2)
idx <- 1
for(i in seq_len(length(term_ids))) {
  for (j in seq(2,length(term_ids))) {
    if (j > i){
      pairwise_terms_1[[idx]] <- term_ids[i]
      pairwise_terms_2[[idx]] <- term_ids[j]
      idx <- idx + 1
    }
  }
}
calc_pairwise_kappa <- function(term_idx, pairwise_terms_1, pairwise_terms_2, go_terms_matrix){
  kappa <- calc_kappa_score(pairwise_terms_1[term_idx], pairwise_terms_2[term_idx], go_terms_matrix)
  return(tibble(
    Term1 = pairwise_terms_1[term_idx],
    Term2 = pairwise_terms_2[term_idx],
    kappa_score = kappa
  ))
}
kappa_scores <- 
  do.call(rbind,
    lapply(seq_len(length(pairwise_terms_1)), calc_pairwise_kappa,
       pairwise_terms_1, pairwise_terms_2, go_terms_matrix))
write_tsv(kappa_scores, path = cmd_line_args$options[['output_file']])
