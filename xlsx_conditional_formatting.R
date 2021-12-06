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
if (!is.null(cmd_line_args$options[['extra_sheets']])) {
  extra_sheets <-
    unlist( strsplit(cmd_line_args$options[['extra_sheets']],
                     ',', fixed = TRUE) )
} else {
  extra_sheets <- NULL
}
yaml_file <- cmd_line_args$options[['yaml']]
guess_max <- cmd_line_args$options[['guess_max']]
debug <- cmd_line_args$options[['debug']]

packages <- c('tidyverse', 'openxlsx', 'yaml')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data
excel_data <- read_tsv(file = cmd_line_args$args[1], na = c("NA"),
                       guess_max = guess_max)

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
    if (is.null(format_info[['type']])) {
      format_type <- "expression"
    } else {
      format_type <- format_info[['type']]
    }
    
    if (is.null(format_info[['rule']])) {
      if (!is.null(format_info[['colname']]) & 
          !is.null(format_info[['operator']]) &
          !is.null(format_info[['text']]) ) {
        # check column exists in data
        if( !(format_info[['colname']] %in% colnames(excel_df)) ) {
          stop(paste0("Column name, ", format_info[['colname']], ", isn't in the data"))
        }
        # work out excel column
        excel_cols <- character(length = ncol(excel_df))
        i <- 0
        j <- 1
        for( colnum in seq_len(ncol(excel_df)) ) {
          excel_cols[colnum] <- paste0(LETTERS[i], LETTERS[j])
          j <- j + 1
          if(j > 26) {
            j <- 1
            i <- i + 1
          }
        }
        excel_colname <- excel_cols[ colnames(excel_df) == format_info[['colname']] ]
        format_rule <- paste0("$", excel_colname, "1", format_info[['operator']], 
                              '"', format_info[['text']], '"')
      } else {
        stop("Don't know how to parse the YAML rules")
      }
    } else {
      format_rule <- format_info[['rule']]
    }
    
    conditionalFormatting(wb, sheet_name, cols=1:ncol(excel_df), 
                          rows=2:(nrow(excel_df)+1), rule=format_rule, 
                          style = createStyle(fontColour = format_info[['fontColour']], 
                                              bgFill = format_info[['bgFill']]),
                          type = format_type)
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
sheet_name <- sub("\\.tsv", "", cmd_line_args$args[1])
if (nchar(sheet_name) > 31) {
  sheet_name <- 'Sheet1'
}

create_and_format_sheet(wb, sheet_name, excel_data, formats)
# create worksheets for extra files
sheet_num <- 2
if (!is.null(extra_sheets)) {
  for (filename in extra_sheets){
    excel_data <- read_tsv(file = filename, na = c("NA"),
                           guess_max = guess_max)
    sheet_name <- sub("\\.tsv", "", filename)
    if (nchar(sheet_name) > 31) {
      sheet_name <- paste0('Sheet', sheet_num)
      sheet_num <- sheet_num + 1
    }
    create_and_format_sheet(wb, sheet_name, excel_data, formats)
  }
}

saveWorkbook(wb, output_file, overwrite = TRUE)
