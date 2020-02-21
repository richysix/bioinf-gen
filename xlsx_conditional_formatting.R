library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='df.xlsx',
              help="Output file name [default %default]" ),
  make_option("--extra_sheets", type="character", default=NULL,
              help="comma-separated list of filenames which are included in the workbook as extra sheets [default %default]" ),
  make_option("--yaml", type="character", default=NULL,
              help="Name of YAML to take formatting rules from [default %default]" ),
  make_option("--guess_max", type="integer", default=1000,
              help="Maximum number of records to use for guessing column types. [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to create an Excel Workbook from a set of tsv files',
  'Conditional formatting of cells can also be included.',
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'xlsx_conditional_formatting.R',
    usage = "Usage: %prog [options] tsv_file",
    description = desc ),
  positional_arguments = 1
)

# cmd_line_args <- list(
#   options = list(output_file = 'df.xlsx',
#                  extra_sheets = '317P2.vep-93.tsv',
#                  yaml = 'var-excel_format.yml',
#                  debug = TRUE),
#   args = c('317P2.vep-93-sig.tsv')
# )

# unpack options
output_file <- cmd_line_args$options[['output_file']]
extra_sheets <-
  unlist( strsplit(cmd_line_args$options[['extra_sheets']],
                   ',', fixed = TRUE) )
yaml_file <- cmd_line_args$options[['yaml']]
guess_max <- cmd_line_args$options[['guess_max']]
debug <- cmd_line_args$options[['debug']]

packages <- c('tidyverse', 'openxlsx', 'yaml')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data
excel_data <- read_tsv(file = cmd_line_args$args[1], na = c("NA"))

# load yaml
if (!is.null(yaml_file)) {
  format_yaml <- read_yaml(file = yaml_file)
  formats <- format_yaml[['formatting']]
} else {
  formats <- NULL
}

create_formatting <- function(sheet_name, formats, excel_df) {
  for(idx in seq_len(length(formats))) {
    format_info <- formats[[idx]]
    if (debug) {
      cat(format_info[['rule']], format_info[['fontColour']], format_info[['bgFill']], "\n")
    }
    conditionalFormatting(wb, sheet_name, cols=1:ncol(excel_df), 
                          rows=1:(nrow(excel_df)+1), rule=format_info[['rule']], 
                          style = createStyle(fontColour = format_info[['fontColour']], 
                                              bgFill = format_info[['bgFill']]) )
  }
}

create_and_format_sheet <- function(workbook, name, excel_df, formats) {
  addWorksheet(workbook, name)
  writeData(workbook, name, excel_df)
  # create formatting statements
  if (!is.null(formats)) {
    create_formatting(name, formats, excel_df)
  }
}

# create workbook
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Calibri")

# create worksheet for main file
create_and_format_sheet(wb, sub("\\.tsv", "", cmd_line_args$args[1]), 
                        excel_data, formats)
# create worksheets for extra files
for (filename in extra_sheets){
  excel_data <- read_tsv(file = filename, na = c("NA"),
                         guess_max = guess_max)
  create_and_format_sheet(wb, sub("\\.tsv", "", filename), 
                          excel_data, formats)
}

saveWorkbook(wb, output_file, overwrite = TRUE)
