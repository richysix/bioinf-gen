#!/usr/bin/env Rscript

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# Based on this script https://github.com/JacobRPrice/JacobRPrice.github.io/blob/master/assets/misc/UpdatedRVersion-Install_Packages.R
# From this blog post https://jacobrprice.github.io/2019/09/19/Installing-multiple-parallel-R-versions.html
# 
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2022. Richard White
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007

library('optparse')

option_list <- list(
  make_option("--get_packages", action = "store_true", type = "logical",
              default = FALSE,
              help = "Get a list of installed packages [default %default]"),
  make_option(
    "--list_packages", action = "store_true", type = "logical",
    default = FALSE,
    help = "List which packages would be installed [default %default]"
  ),
  make_option(
    "--no_check", action = "store_true", type = "logical", default = FALSE,
    help = "Don't check whether packages are already installed [default %default]"
  ),
  make_option(
    c("-d", "--debug"), action = "store_true", type = "logical",
    default = FALSE, help = "print debugging information [default %default]"
  )
)

desc <- paste(
  "Script to install packages from a file of package names",
  "Requires a tab-separated file of package names and which repository to use",
  "To update a previous version of R, first run this script using the previous",
  "version with the --get_packages option.",
  "Then run the script again in the new version to install previously",
  "installed packages. The --list_packages options prints which packages would",
  "be installed.",
  "The first argument should be the packages filename",
  "",
  sep = "\n"
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    prog = 'UpdatedRVersion-Install-Packages.R',
    usage = "Usage: %prog [options] packages_file" ),
  positional_arguments = 1
  )

list_packages <- cmd_line_args$options[['list_packages']]
get_packages <- cmd_line_args$options[['get_packages']]
if (list_packages & get_packages) {
  stop('UpdatedRVersion-Install-Packages.R\n',
       'Only one of --list_packages and --get_packages can be specified')
}
install_packages <- ifelse(list_packages | get_packages, FALSE, TRUE)
no_check <- cmd_line_args$options[['no_check']]

if (cmd_line_args$options[['debug']]) {
  print(cmd_line_args$options)
  cat(
    sprintf('List packages: %s\nGet packages: %s\nInstall packages: %s\n', 
            list_packages, get_packages, install_packages)
  )
}

if (get_packages) {
  if (!require(devtools, quietly = TRUE)) {
    stop("The devtools package must be installed to get the installed packages")
  }
  installed_packages <- 
    installed.packages() |> 
    rownames() |> 
    devtools::package_info() |> 
    as.data.frame()
  write.table(installed_packages[, c("package", "source")],
              file = cmd_line_args$args[1],
              quote = FALSE, sep = "\t",
              row.names = FALSE)
  quit(save = "no")
}
if (install_packages) {
  # set the CRAN mirror
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.ma.imperial.ac.uk/"
  options(repos = r)
  
  print(getOption("repos"))
  # check if devtools is installed and install it always first :-)
  if (!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
  }
  
  # Bioconductor stuff
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install()
}
#######
# now move on to the installation of the most useful items...
###
# define packages we need to install
# open package file
packages <- read.delim(cmd_line_args$args[1])

.cran_packages <- packages$package[ grepl("CRAN", packages$source) ]

.bioc_packages <- packages$package[ grepl("Bioconductor", packages$source) ]

.github_packages <- packages$package[ grepl("Github", packages$source) ]
.github_repos <- packages$source[ grepl("Github", packages$source) ]
.github_repos <- sub("Github \\(", "", .github_repos)
.github_repos <- sub("@[a-z0-9]+\\)", "", .github_repos)

# install packages which are not currently installed.
check_packages_installed <- function(package_list) {
  installed <- logical(length = length(package_list))
  for (i in seq_along(package_list)) {
    if ( length(find.package(package_list[i], quiet = TRUE)) > 0) {
      installed[i] <- TRUE
    }
  }
  return(installed)
}

if (no_check) {
  .inst <- rep(FALSE, length(.cran_packages))
} else {
  .inst <- check_packages_installed(.cran_packages)
}
if (any(!.inst)) {
  if (list_packages) {
    cat("CRAN packages to install:\n", paste(.cran_packages[!.inst], sep = "\n"), "\n")
  } else {
    install.packages(.cran_packages[!.inst])
  }
}

if (no_check) {
  .inst <- rep(FALSE, length(.bioc_packages))
} else {
  .inst <- check_packages_installed(.bioc_packages)
}
if (any(!.inst)) {
  if (list_packages) {
    cat('Bioconductor packages to install:\n', paste(.bioc_packages[!.inst], sep = "\n"), "\n")
  } else {
    BiocManager::install(.bioc_packages[!.inst], update = FALSE)
  }
}

if (no_check) {
  .inst <- rep(FALSE, length(.github_packages))
} else {
  .inst <- check_packages_installed(.github_packages)
}
if (any(!.inst)) {
  if (list_packages) {
    cat('GitHub packages to install:\n', paste(.bioc_packages[!.inst], sep = "\n"), "\n")
  } else {
    devtools::install_github(.github_repos[!.inst])
  }
}
