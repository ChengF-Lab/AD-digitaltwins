#####################################################################
# 0. Load libarary and archR proj
#####################################################################

library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)
addArchRThreads(threads = 32) 
addArchRGenome("hg38")

# Define the input path
atac_dir = "NotAD_out"

##########################################################################################
# 1. Load Previously Prepared ArchR project
##########################################################################################

atac_proj3 <- loadArchRProject(atac_dir, force=TRUE)

print(colnames(getCellColData(atac_proj3)))
table(atac_proj3$NamedClust)

getAvailableMatrices(atac_proj3)


##########################################################################################
# 2. Identifying Marker Peaks with ArchR
##########################################################################################
markersPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj3, 
    useMatrix = "PeakMatrix", 
    groupBy = "NamedClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

## get markers and save it
markerList_all <- getMarkers(markersPeaks)
markerList_all
saveRDS(markerList_all, file = "./marker_peak_all.rds")

markerList_sig <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList_sig
saveRDS(markerList_sig, file = "./marker_peak_sig.rds")

##########################################################################################
# 3. Motif Enrichment
##########################################################################################
atac_proj3 <- addMotifAnnotations(ArchRProj = atac_proj3, motifSet = "cisbp", name = "Motif",force=TRUE)

##########################################################################################
# 3.1 Motif Enrichment in Marker Peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = atac_proj3,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
  )

enrichMotifs  

# Plot a different style heatmap
plot_mat <- plotEnrichHeatmap(enrichMotifs[,column_order], n=5, transpose=FALSE, cutOff=10, returnMatrix=TRUE)

# Extract top 5 enriched genes for each column
top_genes <- lapply(1:ncol(plot_mat), function(i) {
    top_indices <- order(plot_mat[, i], decreasing = TRUE)[1:5]  # Get indices of top 5 values
    return(rownames(plot_mat)[top_indices])  # Return gene names
})

# Create a new matrix with the selected genes in the specific order
new_mat <- do.call(rbind, lapply(1:length(top_genes), function(i) {
    plot_mat[top_genes[[i]], , drop = FALSE]  # Subsetting and maintaining as matrix
}))

# Adjust column and row names if needed
colnames(new_mat) <- colnames(plot_mat)

# Modify row names by removing characters after "_"
rownames(new_mat) <- unlist(top_genes)  # Flattening the list to set row names
rownames(new_mat) <- sub("_.*", "", rownames(new_mat))

saveRDS(new_mat, file = "marker_Topmotif.rds")

library(ComplexHeatmap)
library(circlize)
# Assuming `enrichMotifs` is a matrix or can be converted to one
pdf("Plots/Motifs-Enriched-Marker-Heatmap_Advanced.pdf",width = 6, height = 8)
Heatmap(new_mat, 
        col = colorRampPalette(c("#E6E7E8", "#3A97FF", "#8816A7","#000436"))(100),
        cluster_rows = FALSE,    # Disables clustering of rows
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        name = "Norm.Enrichment -log10(P-adj)[0-Max]",  # Sets the legend title
        row_names_side = "left",    # Puts row names on the left
        border = TRUE,      # Adds a border around the heatmap
        #use_raster = FALSE
        )  
dev.off()

##########################################################################################
# 5. ChromVAR Deviations
##########################################################################################

if("Motif" %ni% names(atac_proj3@peakAnnotation)){
    atac_proj3 <- addMotifAnnotations(ArchRProj = atac_proj3, motifSet = "cisbp", name = "Motif",force = TRUE)
}

## background peaks
atac_proj3 <- addBgdPeaks(atac_proj3)

## compute per-cell deviations accross all of motif annotations
atac_proj3 <- addDeviationsMatrix(
  ArchRProj = atac_proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)

## access these deviations, set plot = TRUE, return a ggplot object
plotVarDev <- getVarDeviations(atac_proj3, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = atac_proj3, addDOC = TRUE)

# Save project
saveArchRProject(atac_proj3)

##########################################################################################
# 8. Peak2GeneLinkage
##########################################################################################

atac_proj3 <- addPeak2GeneLinks(
    ArchRProj = atac_proj3,
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = atac_proj3,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

p2g

markerGenes  <- c(
    "SPIB", "SPI1", "BCL11A", "BCL11B", "SPIC",
    "CTCF","CTCFL","SOX13","SOX9","SOX4",
    "NFIC","NFIX","NFIB","NFIA","ZFX",
    "ZBTB7A","ZNF148","PATZ1",
    "TCF12","NHLH2","ASCL2","ASCL1","TFAP4",
    "JUNB", "JUN","FOSL2","FOS", "JUND"
  )

p <- plotBrowserTrack(
    ArchRProj = atac_proj3, 
    groupBy = "NamedClust", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(atac_proj3)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = atac_proj3, 
    addDOC = FALSE, width = 5, height = 5)


## Plotting a heatmap of peak-to-gene links
p <- plotPeak2GeneHeatmap(ArchRProj = atac_proj3, groupBy = "NamedClust")
p
ggsave("Plots/PlotPeak2GeneHeatmap.pdf")

##########################################################################################
# 9. Identification of Positive TF-Regulators
##########################################################################################
seGroupMotif <- getGroupSE(ArchRProj = atac_proj3, useMatrix = "MotifMatrix", groupBy = "NamedClust")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
    ArchRProj = atac_proj3,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p


library(ggplot2)
library(ggrepel)  # Ensure ggrepel is loaded

# Your existing plot code
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  ) +
  geom_text(data = subset(data.frame(corGSM_MM), GeneScoreMatrix_name %in% markerGenes & TFRegulator == "YES"),
            aes(label = GeneScoreMatrix_name),  # Use your specific label column here
            vjust = -1)  # Adjust text positioning


ggsave("Plots/Pos_TF_corGSM.pdf")


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  ) 

# Filtering and sorting for the text labels
top_labels <- as.data.frame(corGSM_MM) %>%
  filter(TFRegulator == "YES") %>%
  arrange(desc(cor)) %>%
  head(5)

p +  geom_text(data = top_labels,
            aes(label = GeneScoreMatrix_name),  # Use your specific label column here
            vjust = -1)  # Adjust text positioning


ggsave("Plots/Pos_TF_corGSM2.pdf")