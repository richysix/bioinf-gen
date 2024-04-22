#!/usr/bin/env Rscript
desc_string <- "
This script creates a Gantt chart from an Excel spreadsheet
The worksheet must have the following columns
TaskID - Unique id for the task
Task - Label for the task
Category - Label for categories that the task are assigned to
CategoryID - Unique id for the category
Start - Start time of the task
TimeEstimate - Estimated duration of the task

The following columns are optional:
SuperCategory
"

library('optparse')

option_list <- list(
  make_option("--output_file", type = "character", default = 'gantt.pdf',
              help = "Name of output file [default %default]" ),
  make_option("--sheet_name", type = "character", default = NULL, 
              help = paste("Name of sheet to use if it isn't the first sheet",
                           "in the workbook [default %default]" )), 
  make_option("--add_super_category_boxes", type = "logical", default = FALSE, 
              action = "store_true", 
              help = "Draws boxes behind the gantt chart to group supercategories together [default %default]" ), 
  make_option("--width", type = "integer", default = NULL,
              help = "Width of the plot [default %default]" ),
  make_option("--height", type = "integer", default = NULL,
              help = "Height of the plot [default %default]" )
  # make_option("--debug", type = "logical", default = FALSE, 
  #             action = "store_true", 
  #             help = "Turns on debugging statements [default %default]" )
)

# For testing. If running this script interactively the options
# get set to defaults and positional arguments are set to
# whatever is in the arguments vector below
# change this line 
# arguments <- c('data/Amp.sig.tsv', 'data/Oxy.sig.tsv')
# to change the file names
if (any(commandArgs() == "--interactive")) {
  arguments <- c('workplan_data.xlsx')
} else {
  arguments <- commandArgs(trailingOnly = TRUE)
}

# get command line options
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list, prog = 'gantt-from-excel.R', 
    usage = "Usage: %prog [options] input_file", 
    description = desc_string), 
  args = arguments,
  positional_arguments = 1
)

# load packages
packages <- c('tidyverse', 'readxl', 'biovisr', 'miscr')
for (package in packages) {
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data
sheet <- ifelse(is.null(cmd_line_args$options$sheet_name), 1, 
                cmd_line_args$options$sheet_name)
gantt_data <- read_excel(path = cmd_line_args$args[1],
                         sheet = sheet)

required_cols <- c("TaskID", "Task", "CategoryID", "Category", "Start", "TimeEstimate")
if (any(!(required_cols %in% colnames(gantt_data)))) {
  missing_cols <- required_cols[ !(required_cols %in% colnames(gantt_data)) ]
  exit_msg <- paste("Some of the required columns are missing!\n",
                    "Missing colums: ", paste0(missing_cols, collapse = ", "))
  stop(exit_msg)
}

gantt_data <- gantt_data |> 
  mutate(midpoint = Start + TimeEstimate/2,
         End = Start + TimeEstimate)
max_months <- max(gantt_data$End)
max_yr <- ceiling(max_months/12)

# set levels of TaskID and CategoryID
gantt_data <- gantt_data |> 
  mutate(
    TaskID = factor(TaskID, levels = unique(TaskID)),
    CategoryID = factor(CategoryID, levels = unique(CategoryID)),
    label_xpos = case_when(
      End <= max_months/2 ~ End + 1,
      TRUE ~ Start - 1
    ),
    label_hjust = case_when(
      End <= max_months/2 ~ 0,
      TRUE ~ 1
    ))

if (cmd_line_args$options$add_super_category_boxes) {
  # check SuperCategory exists
  if (!("SuperCategory" %in% colnames(gantt_data))) {
    exit_msg <- "Columns must include SuperCategory to add boxes"
    stop(exit_msg)
  }
  boxes <- gantt_data |> 
    group_by(SuperCategory) |> 
    summarise(
      ymax = factor(TaskID[1], levels = rev(levels(TaskID))),
      ymin = factor(TaskID[length(TaskID)], levels = rev(levels(TaskID))),
      End = max(End)
    ) |> 
    mutate(
      ymin = as.integer(ymin) - 0.48,
      ymax = as.integer(ymax) + 0.48,        
      xmin = 0,
      xmax = max_months,
      label_x = case_when(
        End <= max_months/2 ~ xmax - 0.5,
        TRUE ~ xmin + 0.5
      ),
      label_hjust = case_when(
        End <= max_months/2 ~ 1,
        TRUE ~ 0
      )
    )
  
  gantt_plot_by_task <- 
    ggplot(data = gantt_data) +
    geom_tile(aes(x = midpoint, y = TaskID, width = TimeEstimate),
              colour = NA, fill = NA,
              height = 0.95, show.legend = FALSE) +
    geom_rect(data = boxes, show.legend = FALSE, 
              colour = "grey80", fill = NA,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax)) +
    geom_text(data = boxes, colour = "grey30", nudge_y = -0.25, 
              size = 8/.pt, fontface = "bold",
              aes(x = label_x, y = ymax, label = SuperCategory,
                  hjust = label_hjust))
  
} else {
  gantt_plot_by_task <- 
    ggplot(data = gantt_data)
}

yLabels <- gantt_data$Task |> magrittr::set_names(gantt_data$TaskID)
gantt_plot_by_task <- 
  gantt_plot_by_task +
  geom_tile(aes(x = midpoint, y = TaskID,
                fill = CategoryID, width = TimeEstimate),
            height = 0.95, show.legend = FALSE) +
  scale_fill_manual(values = biovisr::cbf_palette(gantt_data$CategoryID),
                    name = "Project Component") + 
  scale_x_continuous(breaks = seq(0, max_months, 12), labels = seq(0, max_yr, 1),
                     name = "Project year") + 
  scale_y_discrete(limits = rev, labels = yLabels) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_blank())

# aggregate to Category
category_summary <- gantt_data |> 
  group_by(CategoryID, Category, SuperCategory) |> 
  summarise(
    End = max(Start + TimeEstimate),
    Start = min(Start),
    .groups = "drop"
  ) |> 
  mutate(
    TimeEstimate = End - Start,
    midpoint = Start + TimeEstimate/2,
    label_xpos = case_when(
      End <= max_months/2 ~ End + 1,
      TRUE ~ Start - 1
    ),
    label_hjust = case_when(
      End <= max_months/2 ~ 0,
      TRUE ~ 1
    )
  )

yLabels <- category_summary$Category |> magrittr::set_names(category_summary$CategoryID)

if (cmd_line_args$options$add_super_category_boxes) {
  boxes <- gantt_data |> 
    group_by(SuperCategory) |> 
    summarise(
      ymax = factor(CategoryID[1], levels = rev(levels(CategoryID))),
      ymin = factor(CategoryID[length(CategoryID)], levels = rev(levels(CategoryID))),
      End = max(End)
    ) |> 
    mutate(
      ymin = as.integer(ymin) - 0.48,
      ymax = as.integer(ymax) + 0.48,        
      xmin = 0,
      xmax = max_months,
      label_x = case_when(
        End <= max_months/2 ~ xmax - 0.5,
        TRUE ~ xmin + 0.5
      ),
      label_hjust = case_when(
        End <= max_months/2 ~ 1,
        TRUE ~ 0
      )
    )
  
  gantt_plot_by_category <- 
    ggplot(data = category_summary) +
    geom_tile(aes(x = midpoint, y = CategoryID, width = TimeEstimate),
              colour = NA, fill = NA,
              height = 0.95, show.legend = FALSE) +
    geom_rect(data = boxes, show.legend = FALSE, 
              colour = "grey80", fill = NA,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax)) +
    geom_text(data = boxes, colour = "grey30", nudge_y = -0.25, 
              size = 8/.pt, fontface = "bold",
              aes(x = label_x, y = ymax, label = SuperCategory,
                  hjust = label_hjust))
  
} else {
  gantt_plot_by_category <- 
    ggplot(data = category_summary)
}

gantt_plot_by_category <- 
  gantt_plot_by_category +
  geom_tile(aes(x = midpoint, y = CategoryID,
                fill = CategoryID, width = TimeEstimate),
            height = 0.95, show.legend = FALSE) +
  scale_fill_manual(values = biovisr::cbf_palette(gantt_data$CategoryID),
                    name = "Project Component") + 
  scale_x_continuous(breaks = seq(0, max_months, 12), labels = seq(0, max_yr, 1),
                     name = "Project year") + 
  scale_y_discrete(limits = rev, labels = yLabels) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())

dimensions <- list(
  "png" = list(width = 1200, height = 600, res = 200),
  "svg" = c(width = 6, height = 3),
  "pdf" = c(width = 6, height = 3),
  "eps" = c(width = 6, height = 3)
)

if (!grepl("pdf$", cmd_line_args$options$output_file)) {
  suffix <- tools::file_ext(cmd_line_args$options$output_file)
  base <- tools::file_path_sans_ext(cmd_line_args$options$output_file)
  output_file <- paste0(base, '-by-task.', suffix)
  params_list <- c(list(NULL), dimensions[[suffix]])
  params_list[[1]] <- list(filename = output_file, plot = gantt_plot_by_task)
  if (!is.null(cmd_line_args$options[["width"]])) {
    params_list$width <- cmd_line_args$options[["width"]]
  }
  if (!is.null(cmd_line_args$options[["height"]])) {
    params_list$height <- cmd_line_args$options[["height"]]
  }
  do.call(miscr::output_plot, params_list)
  
  output_file <- paste0(base, '-by-category.', suffix)
  params_list[[1]] <- list(filename = output_file, plot = gantt_plot_by_category)
  do.call(miscr::output_plot, params_list)
} else {
  output_file <- cmd_line_args$options$output_file
  pdf(file = output_file)
  print(gantt_plot_by_task)
  print(gantt_plot_by_category)
  invisible(dev.off())
}
