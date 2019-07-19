
PackageLoad2 <- function(reqd.packages = NULL){
  
  # Load the required package for using the function: isPackageLoaded() first
  require(R.utils)
  cat("\014")
  
  missing.packages <- reqd.packages[!(reqd.packages %in% as.character(installed.packages()[,'Package']))]
  if(length(missing.packages) > 0) install.packages(missing.packages)
  
  missing.packages <- reqd.packages[!(reqd.packages %in% as.character(installed.packages()[,'Package']))]
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if(length(missing.packages) > 0) lapply(missing.packages, BiocManager::install)
  cat("\014") 
  
  toload.packages <- as.list(reqd.packages[!isPackageLoaded(reqd.packages)])
  lapply(toload.packages, require, character.only = T)
  cat("\014") 
  
}