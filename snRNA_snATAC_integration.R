######################################################
## 1. Load libarary
######################################################

library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(SingleCellExperiment)

set.seed(1)
addArchRThreads(threads = 32) 
addArchRGenome("hg38")

# Define the input path
atac_dir = "snATAC/NotAD_out"
wd_dir = "snATAC/"


######################################################
## 2. Load seRNA and saved snATAC proj
######################################################

# load snATAC data
atac_proj <- loadArchRProject(path = atac_dir) 
# Ensure Impute Weights are up to date
#atac_proj <- addImputeWeights(atac_proj)
colnames(getCellColData(atac_proj))


# load snRNA data
seRNA <- readRDS("Not_AD_sce.rds")
seRNA

colnames(colData(seRNA))
table(colData(seRNA)$celltype2)


######################################################
## Rename the cluster based genescore
######################################################


clustNames <- list(
  "C1" = "Mic",
  "C2" = "Ext", 
  "C3" = "Inh", 
  "C4" = "Oli",
  "C5" = "Ast", 
  "C6" = "OPC"
  )

atac_proj$NamedClust <- clustNames[atac_proj$Clusters] %>% unlist()
table(atac_proj$NamedClust)

saveArchRProject(atac_proj,load = FALSE)


######################################################
## 3. Unconstrained Integration
######################################################

atac_proj2 <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj, 
    force = TRUE,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    sampleCellsATAC = 10000, # Default for both was 10000
    sampleCellsRNA = 10000,
    nGenes = 2000, # Default was 2000
    dimsToUse = 1:40, # Default was 30
    scaleTo = 10000,   # Default was 10000
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    groupRNA = "celltype2", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "predictedCell_Un", # Name of column where cell from scRNA is matched to each cell
    nameGroup = "predictedGroup_Un", # Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore_Un"  # Name of column where prediction score from scRNA
)
atac_proj2
colnames(getCellColData(atac_proj2))
getAvailableMatrices(atac_proj2)

saveArchRProject(atac_proj2)
#saveArchRProject(ArchRProj = atac_proj2, outputDirectory = "snRNA_ATAC_Intergration", load = FALSE)


cM <- as.matrix(confusionMatrix(atac_proj2$NamedClust, atac_proj2$predictedGroup_Un))
cM
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
unique(unique(atac_proj2$predictedGroup_Un))

