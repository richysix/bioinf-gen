#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--key", type="character", default=NULL, action="store",
            help=paste0("Name of the column to use as the key.",
            " Column number if the input has no header [default %default]") ),
  make_option("--value", type="character", default=NULL, action="store",
            help=paste0("Name of the column to use as the value.",
            " Column number if the input has no header [default %default]") ),
  make_option("--fill_value", type="character", default='0', action="store",
            help="Value to use to fill any missing cells [default %default]"),
  make_option("--header", type="logical", default=FALSE, action="store_true",
            help="Does the input data have a header line [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
            help="Turns on debugging statements [default %default]" )
)

desc <- paste('', 'Script to reformat a file from long to wide data',
    'One of the columns in the input data becomes the column names in the output data',
    'The script accepts input on STDIN or from a file. For STDIN use - as the filename',
    sep = "\n")

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'reshape-long_to_wide.R',
    usage = "Usage: %prog [options] inputFile", description = desc ),
  positional_arguments = 1
)
debug <- cmd_line_args$options[['debug']]

library('tidyr')

input_file <- cmd_line_args$args[1]
if (input_file == '-') {
    input_file <- file("stdin")
}
input_data <- read.delim(file = input_file,
                         header = cmd_line_args$options[['header']])

# if key or value columns aren't specified assume the last column is the values 
# and the column before that is the one to spread
check_col <- function(option_name, input_data) {
    option <- cmd_line_args$options[[option_name]]
    if (is.null(option)) {
        if (option_name == 'key') {
            col <- ncol(input_data) - 1
        } else if (option_name == 'value') {
            col <- ncol(input_data)
        } else {
            stop('option_name parameter must be one of "key" or "value"')
        }
    } else {
        # if there is a header assume the key is a column name
        if (cmd_line_args$options[['header']]) {
            # check supplied column name exists
            if (!(option) %in% colnames(input_data)) {
                stop(paste0('The value supplied for option --', option_name,
                            ', "', option, '", is not one of column names\n'))
            } else {
                col <- option
            }
        } else {
            # assume it's a column number
            col <- as.integer(option)
            # check the column exists
            if (col > ncol(input_data)) {
                stop(paste0('The column number supplied for option --',
                           option_name, ', "', option, '", is too large\n'))
            }
        }
    }
    return(col)
}

key <- check_col('key', input_data)
value <- check_col('value', input_data)

if (debug) {
    cat('key =', key,
        '\nvalue =', value, '\n')
}

reshaped_data <- spread(input_data, key, value,
                        fill = cmd_line_args$options[['fill_value']])

write.table(reshaped_data, file = "", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = cmd_line_args$options[['header']])

