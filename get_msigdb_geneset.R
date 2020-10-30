library('optparse')

option_list <- list(
  make_option("--species", type="character", default='Danio rerio',
              help="Species name [default %default]" ),
  make_option("--category", type="character", default=NULL,
              help="MSigDB collection abbreviation, such as H, C1, C2, C3, C4, C5, C6, C7 [default %default]" ),
  make_option("--subcategory", type="character", default=NULL,
              help="MSigDB sub-collection abbreviation, such as CGP or BP [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to download MSigDB genesets',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'get_msigdb_geneset.R',
    usage = "Usage: %prog [options] output_filename",
    description = desc ),
  positional_arguments = 1
)

# cmd_line_args <-
#   list(options = list('output_file' = 'msigdb_genesets.gmt',
#                       debug = TRUE),
#        args = c())

# load packages
packages <- c('msigdbr', 'tidyverse')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# function to turn a geneset data frame into a single row
# in gmt format (columns: name, description, genes)
create_gmt_tibble <- function(df) {
  tibble(
    gs_name = df$gs_name[1],
    gs_desc = sprintf('https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html',
                      df$gs_name[1]),
    members = paste0(df$entrez_gene, collapse = "\t")
  )
}

# retrieve data by species, category and sub-category
m_df <- msigdbr(species = cmd_line_args$options[['species']], 
                category = cmd_line_args$options[['category']], 
                subcategory = cmd_line_args$options[['subcategory']])

# split data by geneset, produce a gmt line per geneset
# and rbind everything back together
gmt <- split(m_df, m_df$gs_id) %>% 
  lapply(., create_gmt_tibble) %>% 
  do.call(rbind, .)

# write out to file
write.table(gmt, file = cmd_line_args$args[1], sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
