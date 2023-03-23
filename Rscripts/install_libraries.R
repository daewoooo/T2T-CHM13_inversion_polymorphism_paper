## Install required libraries ##
################################

## Install Bioconductor packages
BiocManager::install("GenomicRanges")
BiocManager::install("ggtree")
BiocManager::install("regioneR")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
BiocManager::install("ggmsa")

## Install CRAN packages
install.packages("dplyr")
install.packages("parallelDist")
install.packages("ggplot2")
install.packages("phangorn")
install.packages("ape")
install.packages("aplot")

