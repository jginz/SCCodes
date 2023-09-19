########################################################
## Load all the required libraries
## Author: Chaitanya R. Acharya
## Last updated: July 13, 2023
########################################################


## Load all required R libraries
## To install a library, for example, 
##    type: BiocManager::install("Seurat")

### NOTE 1: 
###     Seurat v5 can currently be installed from GitHub
##      remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
##
##      For the following code, please DO NOT run 
##      options(Seurat.object.assay.version = "v5")
##

#########################################################################################
## Check to see if the packages have already been installed and install them if necessary
#########################################################################################
LoadWorkshopLibs = function(){
  # Load required packages
  packages = c("Seurat", "BiocManager", "fgsea", "PCAtools", 
               "dittoSeq", "DoubletFinder", "reticulate", "PCAtools", "plyr", "dplyr", 
               "dittoSeq", "grid", "msigdbr", "fgsea", "ggplot2", "tibble", "HGNChelper", 
               "ggraph", "igraph", "tidyverse", "data.tree", "openxlsx", "escape", "UCell")

  # Install packages that are not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    new_packages <- packages[!installed_packages]
    install.packages(new_packages[new_packages != "BiocManager"])  # Exclude BiocManager from CRAN installation
    if("BiocManager" %in% new_packages) {
      install.packages("BiocManager")
    }
  }

  # Install Bioconductor packages if not installed
  bioconductor_packages <- c("qusage", "fgsea", "msigdbr")
  for(pkg in bioconductor_packages) {
    if (!pkg %in% rownames(installed.packages())) {
      BiocManager::install(pkg)
    }
  }

  # Special case for GitHub package
  if(!"DoubletFinder" %in% rownames(installed.packages())) {
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  }

  # Load packages
  invisible(lapply(packages, library, character.only = TRUE))
}
