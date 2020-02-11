library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='zfs-stages-grcz11.rds',
              help="Working directory [default %default]" ),
  make_option("--sample_names", type="character", default=NULL,
              help="Name of column to use instead of names in count file [default %default]" ),
  make_option("--species", type="character", default='Danio_rerio',
              help="Species name [default %default]" ),
  make_option("--assembly", type="character", default='GRCz11',
              help="Genome assembly [default %default]" ),
  make_option("--ensembl_version", type="integer", default=NULL,
              help="Ensembl version number [default %default]" ),
  make_option("--metadata", type="character", default=NULL,
              help="YAML file of metadata [default %default]" ),
  make_option("--warn_missing_samples", type="logical", action="store_true", default=FALSE,
              help="Warn about missing samples instead of exiting [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'create_deseq_rds.R',
    usage = "Usage: %prog [options] sampleFile countFile\n",
    description = paste("Create SummarizedExpt object for data",
                        "from a count and sample file.",
                        "Metadata can be supplied as options or a YAML file.")
  ),
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list(output_file = 'zfs-stages-grcz11-e99.rds',
#                 metadata = '/lustre/scratch117/maz/team31/resources/zfs-stages/zfs-stages-grcz11-e99.yml',
#                 species = 'Danio_rerio', assembly = "GRCz11", ensembl_version = 99,
#                 sample_names = 'sampleName', warn_missing_samples = FALSE, verbose = FALSE ),
#  args = c('/lustre/scratch117/maz/team31/projects/rnaseq-zfs-stages-grcz11/zfs-rnaseq-grcz11-samples.txt',
#           '/lustre/scratch117/maz/team31/projects/rnaseq-zfs-stages-grcz11/e99/all.tsv')
#)

# load packages
packages <- c('DESeq2', 'yaml', 'rnaseqtools')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# load sample data
samples_file <- cmd_line_args$args[1]
samples <- read.table( samples_file, header=TRUE, row.names=1 )

# check sample_names option
if (!is.null(cmd_line_args$options[["sample_names"]])) {
    # check sample_names option is a column in he samples df
    if ( !(cmd_line_args$options[["sample_names"]] %in% colnames(samples)) ) {
        stop("Supplied samples names column, ",
             cmd_line_args$options[["sample_names"]],
             ", does not appear in the samples file")
    }
}

# order levels
for (col_name in colnames(samples)) {
  samples[[col_name]] <-
    factor(samples[[col_name]],
            levels = unique(samples[[col_name]]))
}

# load count data
data_file <- cmd_line_args$args[2]
#data <- read.delim(data_file, header=TRUE, check.names=FALSE)
data <- load_rnaseq_data(data_file)

# get counts
if (any( grepl('.count$', colnames(data)) )) {
    counts <- data[, grepl('.count$', colnames(data)) &
                      !grepl('normalised.count$', colnames(data)) ]
    colnames(counts) <- gsub('.count$', '', colnames(counts))
    rownames(counts) <- data[, 'GeneID']
} else {
    stop(paste0('Data file, ', data_file, ', does not contain count columns!'))
}

# check for normalised counts
if (any( grepl('normalised.count$', colnames(data)) )) {
    norm_counts <-
        data[, grepl('normalised.count$', colnames(data))]
    colnames(norm_counts) <-
        gsub('.normalised.count$', '', colnames(norm_counts))
    rownames(norm_counts) <- data[, 'GeneID']
} else {
    stop(paste0('Data file, ', data_file, ', does not contain normalised counts!'))
}

# check if any samples are missing from either the samples file or the counts file
in_both <- intersect(rownames(samples), colnames(counts))
samples_only <- setdiff(rownames(samples), colnames(counts))
counts_only <- setdiff(colnames(counts), rownames(samples))
norm_counts_only <- setdiff(colnames(norm_counts), rownames(samples))

# create warning for any missing genes
msg <- list()
if ( length(samples_only) > 0 ) {
    msg[['samples_only']] <- "Some of the samples in the samples file are not present in the counts file.\n"
}

if ( length(counts_only) > 0 ) {
    msg[['counts_only']] <- "Some of the samples in the counts are not present in the samples file.\n"
}

if ( length(norm_counts_only) > 0 ) {
    msg[['norm_counts_only']] <- "Some of the samples in the normalised counts are not present in the samples file.\n"
}

if ( length(msg) > 0 ) {
    if ( cmd_line_args$options[["warn_missing_samples"]] ) {
        warning(msg)
    } else {
      stop(msg)
    }
}

# subset counts/norm_counts to common samples and make sure they are ordered the same
samples <- samples[ in_both, ]
counts <- counts[ , in_both ]
norm_counts <- norm_counts[ , in_both ]

# change sample names if option is specified
if (!is.null(cmd_line_args$options[["sample_names"]])) {
    samples$sampleNamesOrig <- rownames(samples)
    rownames(samples) <- samples[[ cmd_line_args$options[['sample_names']] ]]
    colnames(counts) <- rownames(samples)
    colnames(norm_counts) <- rownames(samples)
}

# row data. chr, start etc for each gene
row_data <- data[, c('GeneID', 'Chr', 'Start', 'End',
                      'Strand', 'Biotype', 'Name', 'Description')]
# Create Summarised Expt object
ExptData <-
    SummarizedExperiment(assays = list(counts = as.matrix(counts),
                                       norm_counts = as.matrix(norm_counts)),
                        rowData = DataFrame(row_data),
                        colData = DataFrame(samples))

if (!is.null(cmd_line_args$options[['metadata']])) {
    metadata <- yaml.load_file(cmd_line_args$options[['metadata']])
    metadata(ExptData) <- metadata
} else if ( any(sapply(c('species', 'assembly', 'ensembl_version'),
                   function(opt){ !is.null(cmd_line_args$options[[opt]]) })) ){
    metadata <- list()
    for(option in c('species', 'assembly', 'ensembl_version')) {
        metadata[[option]] = cmd_line_args$options[[option]]
    }
    metadata(ExptData) <- metadata
}

# write out to rdata file
if (cmd_line_args$options[['verbose']]) {
    cat("Saving RDS file...", "\n")
}
saveRDS(ExptData, compress = "xz",
        file = cmd_line_args$options[['output_file']])
