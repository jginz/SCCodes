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
  packages = c("Seurat","BiocManager", "harmony","devtools", "reticulate", "plyr","dplyr","grid","msigdbr","ggplot2","tibble","HGNChelper","ggraph","igraph","tidyverse","ggtree", "plotly", "data.tree","openxlsx")
  biocpackages = c("DoubletFinder","PCAtools","dittoSeq","fgsea","escape","qusage","UCell", "glmGamPoi")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

  bioc_installed_packages = biocpackages %in% rownames(installed.packages())
  if (any(bioc_installed_packages == FALSE)) {
    BiocManager::install(biocpackages[!bioc_installed_packages])
    }
  
  if(!"DoubletFinder" %in% rownames(installed.packages())) 
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  
  if(!"presto" %in% rownames(installed.packages())) 
    devtools::install_github('immunogenomics/presto')

  
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  invisible(lapply(biocpackages, library, character.only = TRUE))
}
