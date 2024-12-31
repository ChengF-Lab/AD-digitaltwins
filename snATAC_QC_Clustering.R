#####################################
# Step1: creating ArchR proj and QC
#####################################

###  Getting Set Up
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")

Input_dir = "input_fragments"
out_dir = "snATAC_output_dir"


##  InputFiles

inputFiles <- getInputFiles(Input_dir) # fragments_file_dir
inputFiles

###  Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  promoterRegion = c(2000, 100),
  excludeChr = c("chrM", "chrY"),
  nChunk = 5,
  QCDir = "QualityControl"
)

ArrowFiles

### Inferring Doublets

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search, "UMAP"/"LSI".
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "euclidean", verbose = FALSE),
 # outDir = getOutputDirectory(input),
  LSIMethod = 1
)

### Creating an ArchRProject

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = out_dir,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)


## Filter doubletes

proj <- filterDoublets(ArchRProj = proj)


### Saving and Loading an ArchRProject


proj <- saveArchRProject(ArchRProj = proj)


#####################################
# Step2: Clustering
#####################################

### Dimensionality Reduction
proj2 <- addIterativeLSI(
    ArchRProj = proj, 
    iterations = 3,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI",
    sampleCellsPre = 10000,
    sampleCellsFinal = 30000,
    clusterParams = list(resolution = c(0.5),sampleCells = 10000, n.start= 10),
    varFeatures = 25000,
    force = TRUE,                        
    seed =1
    )

### call clusters
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE,
    sampleCells = 10000,
    maxClusters =6,
    seed = 1
    )
### Visualizing in a 2D UMAP Embedding
proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI",
    name = "UMAP", 
    nNeighbors = 50, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE,
    seed =1
    )    
##  add imputation weights using MAGIC
proj2 <- addImputeWeights(proj2)
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", 
       name = "Clusters", embedding = "UMAP")
plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj2, addDOC = TRUE, width = 5, height = 5)


proj2 <- saveArchRProject(ArchRProj = proj2)

#####################################
# Step3: Marker Peaks
#####################################

Identifying Marker Genes

markersGS <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
saveRDS(markerList, file = marker_output_file)


### (7-1) selected Marker Genes
markerGenes  <- c(
                  
    'DOCK8','APBB1IP','CD74','CSF1R', 'ARHGAP15', # Mic
    'FLT1', 'EPAS1', 'ABCB1','CLDN5', # Endo
    'EBF1','LAMA2', # VLMC
    'MOG','MOBP','CLDN11','PLP1', # Oli
    'AQP4', 'GFAP', 'ATP1B2', 'GJA1', 'GPC5','RYR3', # Ast
    'CSPG4', 'PDGFRA','OLIG1','SOX10','VCAN', # OPC
    'ADARB2','NXPH1','GAD1','GAD2', # Inh
    'CBLN2', 'TSHZ2', 'SATB2' # Ext
  )

###  Heatmap plot
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  #nLabel = 1,
  #nPrint = 1,
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)


### Session Information

Sys.Date()
sessionInfo()    