#!/usr/bin/env Rscript

#####################################################################
# Call peaks
#####################################################################

library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library("BSgenome.Hsapiens.UCSC.hg38")

set.seed(1)
addArchRThreads(threads = 32) 
addArchRGenome("hg38")

# Define the input path
atac_dir = "NotAD_out"
wd_dir = "snATAC/"


##########################################################################################
# Load Previously Prepared ArchR project
##########################################################################################

atac_proj2 <- loadArchRProject(atac_dir, force=TRUE)

##########################################################################################
# Call Peaks
##########################################################################################

# Create Group Coverage Files that can be used for downstream analysis
atac_proj2 <- addGroupCoverages(
    ArchRProj = atac_proj2, 
    groupBy = "NamedClust",
    minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
    force=TRUE)

# Find Path to Macs2 binary
pathToMacs2 <- findMacs2()

# Call Reproducible Peaks
atac_proj2 <- addReproduciblePeakSet(
    ArchRProj = atac_proj2, 
    groupBy = "NamedClust", 
    maxPeaks = 550000,
    pathToMacs2 = pathToMacs2,
    force = TRUE
)

getGeneAnnotation(atac_proj2)
getGeneAnnotation(atac_proj2)$genes

# Add Peak Matrix
atac_proj3 <- addPeakMatrix(atac_proj2)

# view archR proj infor
getPeakSet(atac_proj3)
getAvailableMatrices(atac_proj3)

# Save project
saveArchRProject(atac_proj3,load = FALSE)
