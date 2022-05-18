#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='kappa_scores.tsv',
              help="Output file name [default %default]" ),
  make_option("--pval_colname", type="character", default='p value',
              help="Name of the column which indicates which genes are significant [default %default]" ),
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
    usage = "Usage: %prog [options] topgo_sig_genes_file",
    description = desc ),
  positional_arguments = 1
)

# load packages
packages <- c('tidyverse', 'rlang')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# topgo sig.genes.tsv
go_sig_genes <- read_tsv(cmd_line_args$args[1])

# check sig col
pval_colname <- rlang::sym(cmd_line_args$options[['pval_colname']])
if (!any(colnames(go_sig_genes) == pval_colname)) {
  stop("Significance column name provided by --pval_colname doesn't match any of the column names")
}

# The pvalue column is used to set TRUE/FALSE on gene to go term
# This is a matrix of genes to GO terms
go_terms_matrix <- 
  go_sig_genes %>% 
  mutate(., pvalue = case_when(
    is.na(!!pval_colname) ~ FALSE, 
    !!pval_colname == 0 ~ FALSE,
    !!pval_colname == 1 ~ TRUE)
  ) %>% 
  select(., Gene, GO.ID, pvalue) %>% 
  pivot_wider(., names_from = GO.ID, values_from = pvalue,
              values_fill = list(pvalue = FALSE)) 

calc_kappa_score <- function(first_col, second_col, go_terms_matrix){
  num_genes <- nrow(go_terms_matrix)
  InBoth <- sum(go_terms_matrix[[first_col]] & go_terms_matrix[[second_col]])
  InFirstNotSecond <- sum(go_terms_matrix[[first_col]] & !go_terms_matrix[[second_col]])
  InSecondNotFirst <- sum(!go_terms_matrix[[first_col]] & go_terms_matrix[[second_col]])
  InNeither <- sum(!go_terms_matrix[[first_col]] & !go_terms_matrix[[second_col]])
  po <- (InBoth + InNeither) / num_genes
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
write_tsv(kappa_scores, file = cmd_line_args$options[['output_file']])

# cmd_line_args for testing
# cmd_line_args <-
#   list(options = list('output_file' = 'kappa-scores/zmp_ph71.kappa_scores.tsv',
#                       'gene_id_colname' = "e92 Ensembl Gene ID",
#                       'pval_colname' = "Sig?",
#                       debug = TRUE),
#        args = c('hom_vs_het_wt/sig.tsv', 'hom_vs_het_wt.overlap/GO.sig.genes.tsv'))

