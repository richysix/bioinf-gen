library('optparse')

option_list <- list(
  make_option("--sig_genes", type="character", action="store", default=NULL,
              help="File of significant genes (in first column) to subset to [default %default]")
)

desc <- paste('', 'Script to convert a file mapping gene ids to GO Terms',
              '', sep = "\n")

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'create_geneset_file.R',
    usage = "Usage: %prog [options] input_file",
    description = desc),
  positional_arguments = 1
)

library('tidyverse')

# read in data
input_data <- read_tsv('danio_rerio_e98_go.txt', col_names = FALSE)

# subset genes to ones from sig genes file
if (!is.null(cmd_line_args$options[['sig_genes']])) {
  # load sig_genes
  sig_genes <- read_tsv(cmd_line_args$options[['sig_genes']],
                        col_names = FALSE)
  # function to return the rows for a gene
  go_terms_for <- function(gene, input_data) {
    filter(input_data, X1 == gene)
  }
  # run go_terms_for over all sig_genes
  sig_genes_subset <- lapply(sig_genes[[1]], go_terms_for, input_data)
  # turn the list into a data frame
  input_data <- do.call(rbind, sig_genes_subset)
}

# function to return a list of gene ids for a GO term
genes_for <- function(go_term, input_data) {
  filter(input_data, X2 == go_term) %>% 
    select(., X1) %>% pull()
}

# run genes_for function over all unique GO terms
go2genes <- lapply(unique(input_data[[2]]), genes_for, input_data)
names(go2genes) <- unique(input_data[[2]])

# loop through list and output
cat(file = 'danio_rerio_e98.gmt',
    append = FALSE)
for (index in seq_len(length(go2genes))) {
  cat(names(go2genes)[index], "\t",
      paste(go2genes[[index]], collapse = "\t"),
      "\n",
      file = 'danio_rerio_e98.gmt',
      append = TRUE)
}
