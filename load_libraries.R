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
  # load libraries
  packages = c("Seurat","DoubletFinder","PCAtools","plyr","dplyr","dittoSeq","grid","msigdbr","fgsea","ggplot2","tibble","HGNChelper","ggraph","igraph","tidyverse", "data.tree","openxlsx","escape","qusage","UCell")

  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  if(!"DoubletFinder" %in% rownames(installed.packages())) 
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
}
