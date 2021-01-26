#!/usr/bin/env Rscript

library('optparse')

#### Please make sure that Cytoscape >v3.6+ is started and the Cytoscape Apps 'yFiles Layout Algorithms' and 'ClueGO v2.5.2' are installed before running this script! ####

option_list <- list(
  make_option("--output_image_file", type="character", default='go_network.svg', 
              help="Output file name (svg, png or pdf) [default %default]" ), 
  make_option("--output_scale", type="numeric", default=NULL, 
              help="Scale factor of network for output image [default Fit to screen]" ), 
  make_option("--output_basename", type="character", default='cluego', 
              help="Basename for output files [default %default]" ), 
  make_option("--organism", type="character", default='Danio rerio', 
              help="Version of ClueGO to use [default %default]" ), 
  make_option("--analysis_name", type="character", default='cluego', 
              help="Name of the analysis [default %default]" ), 
  make_option("--destroy_network", action="store_true", type="logical", default=FALSE, 
              help="Destroy network at the end [default %default]" ), 
  make_option("--min_go_tree_level", type="integer", default=3, 
              help="Lowest (Least specific) level of GO terms to include [default %default]" ), 
  make_option("--max_go_tree_level", type="integer", default=8, 
              help="Highest (Most specific) level of GO terms to include [default %default]" ), 
  make_option("--pvalue_threshold", type="numeric", default=0.05, 
              help="Adjusted pvalue threshold [default %default]" ), 
  make_option("--kappa_score_level", type="numeric", default=0.4, 
              help="Kappa score level for merging groups [default %default]" ), 
  make_option("--get_ontologies", type="logical", default=FALSE, action="store_true",
              help="Return possible ontologies for species and exit [default %default]" ),
  make_option("--verbose", type="logical", default=FALSE, action="store_true", 
              help="Turns on verbose output [default %default]" ), 
  make_option("--debug", type="logical", default=FALSE, action="store_true", 
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to run GO enrichment using ClueGO', 
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'run_cluego.R', 
    usage = "Usage: %prog [options] input_files", 
    description = desc ), 
  positional_arguments = c(1, 8)
)

# cmd_line_args <-
#   list(options = list('output_image_file' = 'go_network.svg',
#                       "output_basename" = 'cluego',
#                       'organism' = 'Danio rerio',
#                       "analysis_name" = 'cluego',
#                       "destroy_network" = FALSE,
#                       "min_go_tree_level" = 3,
#                       "max_go_tree_level" = 8,
#                       "pvalue_threshold" = 0.05,
#                       "kappa_score_level" = 0.4,
#                       verbose = TRUE,
#                       debug = TRUE),
#        args = c('sig-genes.tsv'))
 
# unpack options
debug <- cmd_line_args$options[['debug']]
verbose <- cmd_line_args$options[['verbose']]
output_basename <- cmd_line_args$options[['output_basename']]

# load packages
packages <- c('tidyverse', 'xml2', 'RJSONIO', 'httr', 'biovisr')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# Helper function to transform output into a data frame
text_to_data_frame <- function(table.text) {
  table <- NULL
  rows <- unlist(strsplit(table.text, split="\n"))
  header <- t(unlist(strsplit(rows[1], split="\t")))
  for(i in 2:length(rows)) {
    if(is.null(table)) {
      table <- t(unlist(strsplit(rows[i], split="\t")))
    } else {
      table <- rbind(table, t(unlist(strsplit(rows[i], split="\t"))))
    }
  }
  table <- as.data.frame(table)
  names(table) <- header
  return(table)
}

port_number <- 1234
host_address <- "localhost"

# define base urls
cytoscape_base_url = paste("http://", host_address, ":", toString(port_number), "/v1", sep="")
cluego_base_url = paste(cytoscape_base_url, "apps", "cluego", "cluego-manager", sep="/")

analysis_name <- cmd_line_args$options[['analysis_name']]
if (verbose) {
  print(paste("Analysis: ", analysis_name, sep=""))
  print(paste("Cytoscape Base URL: ", cytoscape_base_url, sep=""))
  print(paste("ClueGO Base URL: ", cluego_base_url, sep=""))
}

#### 0.0 Start up ClueGO in case it is not running yet
response <- POST(url=paste(cytoscape_base_url, "apps", "cluego", "start-up-cluego", sep="/"), encode = "json")
# wait 2 seconds to make sure ClueGO is started
if(http_status(response)$category == "Server error") {
  if (verbose) {
    print("wait 2 secs")
  }
  Sys.sleep(2)
}

#### 1.0 Select the ClueGO Organism to analyze ####
organism_name = cmd_line_args$options[['organism']] # (run "1.1 Get all ClueGO organisms" to get all options)
if (verbose) {
  print(paste("1.0 Select the ClueGO Organism to analyze: ", organism_name, sep=""))
}
response <- PUT(url=paste(cluego_base_url, "organisms", "set-organism", 
                          URLencode(organism_name), sep="/"), encode = "json")
stop_for_status(response, "Set species")

# print all available ontologies and exit if --get_ontologies set
if (cmd_line_args$options[['get_ontologies']]) {
  response <- GET(url=paste(cluego_base_url, "ontologies", "get-ontology-info", sep="/"), encode = "json")
  stop_for_status(response, "Get Ontologies")
  content_list <- content(response, encode = "json")
  ontologies <- sapply(names(content_list), 
                       function(x, content_list){ sprintf("Index = %s, %s", x, content_list[[x]]) },
                       content_list)
  names(ontologies) <- NULL
  cat(ontologies, sep = "\n")
  quit(save = "no")
}

## [optional functions and settings, un comment and modify if needed]
# 1.1 Get all ClueGO organisms
# response <- GET(paste(cluego_base_url, "organisms", "get-all-installed-organisms", sep="/"))
# stop_for_status(response, "Get organisms")
# print(content(response))
#
# 1.2 Get all info for installed organisms
# response <- GET(paste(cluego_base_url, "organisms", "get-all-organism-info", sep="/"))
# stop_for_status(response, "Get all organism info")
# print(content(response))

#### 2.0 Upload IDs for a specific Cluster ####
if (verbose) {
  print(paste("2.0 Upload IDs for Clusters", sep=""))
}

# [optional functions and settings, un comment and modify if needed]
# 2.1 Select the ClueGO ID type
id_type_name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options) "# Automatic #"
response <- PUT(url=paste(cluego_base_url, "ids", "set-id-type", id_type_name, sep="/"), encode = "json")
stop_for_status(response, "Set ID type")

# # 2.2 Refresh ClueGO source files
# POST(url=paste(cluego_base_url, "ids", "refresh-cluego-id-files", sep="/"), encode = "json")
# 
# # 2.3 Get ClueGO ID types
# response <- GET(paste(cluego_base_url, "ids", "get-all-installed-id-types", sep="/"))
# stop_for_status(response, "Get all ID types")
# print(content(response))

# load gene lists
# 2.4 Set the number of Clusters
number_clusters <- length(cmd_line_args$args)
response <- PUT(url=paste(cluego_base_url, "cluster", "max-input-panel", number_clusters, sep="/"),  encode = "json")
stop_for_status(response, "Set Number of Clusters")

if (number_clusters == 2) {
  colour_palette <- c("#A83232", "#323CA8")
} else {
  colour_palette <- biovisr::cbf_palette(number_clusters)
}

# 2.5 Set analysis properties for a Cluster, repeat for each input cluster
# Possible node shapes are "Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V"
cluster_node_shapes <- c("Round Rectangle", "Octagon", "Parallelogram", 
                         "Triangle", "V", "Rectangle", "Diamond", "Hexagon")
for (cluster_index in seq_len(length(cmd_line_args$args))) {
  node_shape = cluster_node_shapes[cluster_index]
  cluster_color = colour_palette[cluster_index]
  min_number_of_genes_per_term = 3 # defaults. could have as an option
  min_percentage_of_genes_mapped = 4 # defaults. could have as an option
  no_restrictions = FALSE # TRUE for no restricions in number and percentage per term
  response <- PUT(url=paste(cluego_base_url, "cluster", "set-analysis-properties", as.character(cluster_index), 
                            URLencode(node_shape, reserved = TRUE), URLencode(cluster_color, reserved = TRUE), 
                            min_number_of_genes_per_term, min_percentage_of_genes_mapped, 
                            no_restrictions, sep="/"), encode = "json")
  stop_for_status(response, "Set Cluster properties")
}

# go through again and load gene lists
identifiers <- list(
  'Homo Sapiens' = 'ENSG[0-9]{11}',
  'Mus Musculus' = 'ENSMUSG[0-9]{11}',
  'Danio rerio' = 'ENSDARG[0-9]{11}'
)
for (cluster_index in seq_len(length(cmd_line_args$args))) {
  file <- cmd_line_args$args[cluster_index]
  gene_info <- read_tsv(file)
  # check gene ids match
  if (is.null(identifiers[[organism_name]])) {
    warning("No identifier regex for species, ", organism_name)
  }else {
    if (!all(grepl(identifiers[[organism_name]], gene_info$Gene))) {
      stop(paste0("Gene ids do not match expected format: species = ", organism_name, 
                  ", id format = ", identifiers[[organism_name]], "\n",
                  "ids = ", paste0(head(gene_info$Gene), collapse = ", ")))
    }
  }
  gene_list <- toJSON(gene_info$Gene)
  response <- PUT(url=paste(cluego_base_url, "cluster", "upload-ids-list", URLencode(as.character(cluster_index)), sep="/"), 
                  body=gene_list, encode = "json", content_type_json())
  stop_for_status(response, "Upload gene lists")
}

# 2.6 Select visual style
# Possibles (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
# start with visual.style = "ShowGroupDifference"
visual_style = "ShowGroupDifference"
response <- PUT(url=paste(cluego_base_url, "cluster", "select-visual-style", visual_style, sep="/"), encode = "json")
stop_for_status(response, "Set visual style")

####  3.0 Select Ontologies
if (verbose) {
  print(paste("3.0 Select Ontologies", sep=""))
}
# Possible node shapes are "Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V"
selected_ontologies <- toJSON(c("1;Ellipse", "2;Hexagon", "4;Rectangle", "5;Diamond")) # (run "3.1 Get all available Ontologies" to get all options)
response <- PUT(url=paste(cluego_base_url, "ontologies", "set-ontologies", sep="/"), 
                body=selected_ontologies, encode = "json", content_type_json())
stop_for_status(response, "Select ontologies and set shape")

## [optional functions and settings, un comment and modify if needed]
# 3.1 Get all available Ontologies
# response <- GET(paste(cluego_base_url, "ontologies", "get-ontology-info", sep="/"))
# stop_for_status(response, "Get all ontologies")
# print(content(response))
# 
# # 3.2 Select Evidence Codes
# evidence.codes <- toJSON(c("All")) # (run "3.3 Get all available Evidence Codes" to get all options)
# response <- PUT(url=paste(cluego_base_url, "ontologies", "set-evidence-codes", sep="/"), body=evidence.codes, encode = "json", content_type_json())
# stop_for_status(response, "Set evidence codes")
#
# # 3.3 Get all available Evidence Codes
# response <- GET(paste(cluego_base_url, "ontologies", "get-evidence-code-info", sep="/"))
# stop_for_status(response, "Get all evidence codes")
# print(content(response))
# 
# 3.4 Set min, max GO tree level
min_level = cmd_line_args$options[['min_go_tree_level']]
max_level = cmd_line_args$options[['max_go_tree_level']]
all_levels = FALSE
response <- PUT(url=paste(cluego_base_url, "ontologies", "set-min-max-levels", 
                          min_level, max_level, all_levels, sep="/"), encode = "json")
stop_for_status(response, "Set GO levels")

# 3.5 Use GO term significance cutoff
p_value_cutoff = cmd_line_args$options[['pvalue_threshold']]
use_significance_cutoff = TRUE
response <- PUT(url=paste(cluego_base_url, "ontologies", use_significance_cutoff, p_value_cutoff, sep="/"), encode = "json")
stop_for_status(response, "Set pvalue threshold")

# # 3.6 Use GO term fusion
use_go_term_fusion = FALSE
response <- PUT(url=paste(cluego_base_url, "ontologies", use_go_term_fusion, sep="/"), encode = "json")
stop_for_status(response, "Set GO term fusion")

# 3.7 Set statistical parameters
enrichment_type = "Enrichment (Right-sided hypergeometric test)" # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
multiple_testing_correction = "Bonferroni step down" # ("Bonferroni", "Bonferroni step down", "Benjamini-Hochberg", "None")
use_mid_pvalues = FALSE
use_doubling = FALSE
useCustomReferenceSet = FALSE
fullPathToReferenceSet = "none"
response <- PUT(url=paste(cluego_base_url, "stats", URLencode(enrichment_type, reserved = TRUE), 
                          URLencode(multiple_testing_correction, reserved = TRUE), 
                          use_mid_pvalues, use_doubling, 
                          useCustomReferenceSet, fullPathToReferenceSet, sep="/"), 
                encode = "json")
stop_for_status(response, "Set stats parameters")

# 3.8 Set the Kappa Score level
kappa_score = cmd_line_args$options[['kappa_score_level']]
response <- PUT(url=paste(cluego_base_url, "ontologies", "set-kappa-score-level", kappa_score, sep="/"), encode = "json")
stop_for_status(response, "Set kappa score level")

# 3.9 Set grouping parameters
do_grouping = TRUE
coloring_type = "Fix" # ("Random", "Fix")
group_leading_term = "Highest Significance" # ("Highest Significance", "#Genes / Term", "%Genes / Term", "%Genes / Term vs Cluster")
grouping_type = "Kappa Score" # ("Kappa Score", "Tree")
init_group_size = 1
perc_groups_for_merge = 50
perc_terms_for_merge = 50
response <- PUT(url=paste(cluego_base_url, "grouping", do_grouping, coloring_type, 
                          URLencode(group_leading_term, reserved = TRUE), 
                          URLencode(grouping_type, reserved = TRUE), 
                          init_group_size, perc_groups_for_merge, perc_terms_for_merge, sep="/"), encode = "json")
stop_for_status(response, "Set grouping prameters")

#### 4.0 Run ClueGO Analysis ####
if (verbose) {
  print(paste("4.0 Run ClueGO Analysis", sep=""))
}
# Run the analysis and save log file
analysis_option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
response <- GET(paste(cluego_base_url, URLencode(analysis_name), URLencode(analysis_option), sep="/"))
if (response$status_code == 404) {
  print(content(response, encode = "json")[["message"]])
}
stop_for_status(response, "Run analysis")
log_file_name = paste(output_basename, "analysis-log.txt", sep = "-")
writeLines(content(response, encoding="UTF-8"), log_file_name)
# print(content(response, encode = "text"))

# 4.1 Get network id (SUID) (CyRest function from Cytoscape)
response <- GET(paste(cytoscape_base_url, "networks", "currentNetwork", sep="/"))
stop_for_status(response, "Get network id")
current_network_suid <- content(response, encode = "json")$data$networkSUID
# print(current_network_suid)

# hide the control and table panels
response <- GET(paste(cytoscape_base_url, "ui", "panels", sep="/"))
stop_for_status(response, "Get panel info")
panel_info <- content(response, encode = "text", encoding = "UTF-8")
# change panels
for (idx in seq_len(length(panel_info))) {
  if (panel_info[[idx]]$name == "SOUTH" | panel_info[[idx]]$name == "WEST") {
    panel_info[[idx]]$state <- "HIDE"
  }
}
# set panel info
response <- PUT(url=paste(cytoscape_base_url, "ui", "panels", sep="/"), 
                body=toJSON(panel_info), encode = "json", content_type_json())
stop_for_status(response, "Set panel info")

# set scale factor
# fit network to current window
response <- GET(paste(cytoscape_base_url, "apply", "fit", current_network_suid, sep="/"))
stop_for_status(response, "Fit network to window")
# get view id
response <- GET(url = paste(cytoscape_base_url, "networks", "views", "currentNetworkView", sep = "/"))
stop_for_status(response, "Get View id")
view_id <- content(response, encode = "json")$data$networkViewSUID

if (is.null(cmd_line_args$options[['output_scale']])) {
  # make 10% smaller to try and fit text in
  # get current scale factor
  response <- GET(url = paste(cytoscape_base_url, "networks", current_network_suid, "views", 
                              view_id, "network", sep = "/"))
  stop_for_status(response, "Get Visual Properties for View")
  current_vis_prop <- content(response, encode = "json")
  current_scale_factor <- unlist(
    lapply(current_vis_prop, function(x){ if (x$visualProperty == "NETWORK_SCALE_FACTOR"){ return(x$value) } }))
  if (verbose) {
    print(paste("Scale factor to fit network:", current_scale_factor))
  }
  # set scale factor to 90% of current
  vis_prop <- list( list("visualProperty" = "NETWORK_SCALE_FACTOR", 
                         "value" = round(0.9 * current_scale_factor, digits = 2)) )
} else {
  # set scale factor to output_scale option
  vis_prop <- list( list("visualProperty" = "NETWORK_SCALE_FACTOR", 
                         "value" = cmd_line_args$options[['output_scale']]) )
}
response <- PUT(url=paste(cytoscape_base_url, "networks", current_network_suid, "views", 
                          view_id, "network", sep = "/"), bypass = "false",
                body = toJSON(vis_prop), encode = "json", content_type_json())
stop_for_status(response, "Change scale factor")

if (verbose) {
  print(paste("Save results", sep=""))
}

# Get network graphics (CyRest function from Cytoscape)
Sys.sleep(1)
image_type = sub("^.*\\.", "", cmd_line_args$options[['output_image_file']]) # svg, png, pdf
if (!(image_type %in% c('svg', 'png', 'pdf'))) {
  stop(paste("Output image must be either 'svg', 'png' or 'pdf' -", 
             cmd_line_args$options[['output_image_file']]))
}
response <- GET(paste(cytoscape_base_url, "networks", current_network_suid, "views", paste("first.", image_type, sep=""), sep="/"))
stop_for_status(response, "First view")
writeBin(content(response, as = "raw"), cmd_line_args$options[['output_image_file']])

# 4.2 Get ClueGO result table
response <- GET(paste(cluego_base_url, "analysis-results", "get-cluego-table", current_network_suid, sep="/"))
stop_for_status(response, "Cluego results table")
result_table_text <- content(response, encode = "text", encoding = "UTF-8")
table_file_name = paste0(output_basename, "-result-table.txt")
write.table(text_to_data_frame(result_table_text), file=table_file_name, row.names=FALSE, quote = FALSE, na="", col.names=TRUE, sep="\t")

# # 4.3 Get ClueGO genes and main functions
# number_of_functions_to_add = 3
# response <- GET(paste(cluego_base_url, "analysis-results", "get-main-functions", current_network_suid, number_of_functions_to_add, sep="/"))
# stop_for_status(response, "Get genes and functions")
# result_table_text <- content(response, encode = "text", encoding = "UTF-8")
# table_file_name = paste0(output_basename, "-genes-with-main-functions.txt")
# write.table(text_to_data_frame(result_table_text), file=table_file_name, row.names=FALSE, quote = FALSE, na="", col.names=TRUE, sep="\t")

# # 4.4 Get genes result table
# response <- GET(paste(cluego_base_url, "analysis-results", "get-gene-table", current_network_suid, sep="/"))
# stop_for_status(response, "Get gene table")
# result_table_text <- content(response, encode = "text", encoding = "UTF-8")
# table_file_name = paste0(output_basename, "-gene-table.txt")
# write.table(text_to_data_frame(result_table_text), file=table_file_name, row.names=FALSE, quote = FALSE, na="", col.names=TRUE, sep="\t")

# 4.5 Get Kappascore Matrix
response <- GET(paste(cluego_base_url, "analysis-results", "get-kappascore-matrix", current_network_suid, sep="/"))
stop_for_status(response, "Get kappa score matrix")
result_table_text <- content(response, encode = "text", encoding = "UTF-8")
table_file_name = paste0(output_basename, "-kappascore-matrix.txt")
write.table(text_to_data_frame(result_table_text), file=table_file_name, row.names=FALSE, quote = FALSE, na="", col.names=TRUE, sep="\t")

# 4.6 Get binary Gene-Term Matrix
response <- GET(paste(cluego_base_url, "analysis-results", "get-binary-gene-term-matrix", current_network_suid, sep="/"))
stop_for_status(response, "Get binary Gene-Term Matrix")
result_table_text <- content(response, encode = "text", encoding = "UTF-8")
table_file_name = paste0(output_basename, "-binary-gene-term-matrix.txt")
write.table(text_to_data_frame(result_table_text), file=table_file_name, row.names=FALSE, quote = FALSE, na="", col.names=TRUE, sep="\t")

# # 4.7 ClueGO Result Chart
# # Get result charts for both cluster as pie chart
# chart.type = "PieChart" # ("PieChart", "BarChart")
# image.type = "svg" # ("svg", "png", "pdf")
# response <- GET(paste(cluego_base_url, "analysis-results", "get-cluego-result-chart", current_network_suid, cluster1, chart.type, image.type, sep="/"))
# stop_for_status(response, "")
# image.file.name = paste(output.folder, paste("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", image.type, sep=""), sep="/")
# writeBin(content(response, encode = "raw"), image.file.name)
# if(example.selection == "ClueGO Rest Example for two gene lists") {
#   response <- GET(paste(cluego_base_url, "analysis-results", "get-cluego-result-chart", current_network_suid, cluster2, chart.type, image.type, sep="/"))
#   stop_for_status(response, "")
#   image.file.name = paste(output.folder, paste("ClueGO-", chart.type, "-For-Cluster", cluster2, ".", image.type, sep=""), sep="/")
#   writeBin(content(response, encode = "raw"), image.file.name)
# }
# 
## [optional functions and settings, un comment and modify if needed]
# # 4.8 Remove ClueGO analysis result
# print(paste("Remove ClueGO Network", sep=""))
# # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
# response <- DELETE(paste(cluego_base_url, "remove-cluego-analysis-result", current_network_suid, sep="/"))
# stop_for_status(response, "")
#
# 4.9 Remove all ClueGO analysis results
# print(paste("Remove all ClueGO Networks", sep=""))
# response <- DELETE(paste(cluego_base_url, "remove-all-cluego-analysis-results", sep="/"))
# stop_for_status(response, "")

# Output network with cluster differences
if (number_clusters > 1) {
  visual_style = "ShowClusterDifference"
  response <- PUT(url=paste(cluego_base_url, "cluster", "select-visual-style", visual_style, sep="/"), encode = "json")
  stop_for_status(response, "Set visual style")
  Sys.sleep(5)
  
  image.type = sub("^.*\\.", "", cmd_line_args$options[['output_image_file']]) # svg, png, pdf
  if (!(image.type %in% c('svg', 'png', 'pdf'))) {
    stop(paste("Output image must be either 'svg', 'png' or 'pdf' -", 
               cmd_line_args$options[['output_image_file']]))
  }
  response <- GET(paste(cytoscape_base_url, "networks", current_network_suid, "views", paste("first.", image.type, sep=""), sep="/"))
  stop_for_status(response, "Get first view")
  cluster_image_file <- sub("(^.*)\\.([a-z])", "\\1-cluster.\\2", 
                           cmd_line_args$options[['output_image_file']])
  writeBin(content(response, as = "raw"), cluster_image_file)
}

if (verbose) {
  print(paste0("Analysis ", analysis_name, " done"))
}

# # names for all visual styles
# response <- GET(paste(cytoscape_base_url, "apply", "styles", sep = "/"))
# stop_for_status(response, "Get visual styles info")
# visual_styles <- content(response, encode = "text", encoding = "UTF-8")
# style_name <- "ClueGOVisualStyleForGroups_0"
# 
# # mappings
# response <- GET(paste(cytoscape_base_url, "styles", style_name, "mappings", sep = "/"))
# stop_for_status(response, "")
# style_mappings <- content(response, encode = "text", encoding = "UTF-8")

if (cmd_line_args$options[['destroy_network']]) {
  # make sure everything has finished
  Sys.sleep(5)
  # 4.8 Remove ClueGO analysis result
  if (verbose) {
    print(paste("Remove ClueGO Network", sep=""))
  }
  # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
  response <- DELETE(paste(cluego_base_url, "remove-cluego-analysis-result", current_network_suid, sep="/"))
  stop_for_status(response, "Delete network")
}
