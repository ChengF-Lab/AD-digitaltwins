library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(ExperimentHub)
#library("zellkonverter")
library(scater)
library(data.table)
library(SingleCellExperiment)
library(UpSetR)
library(cowplot)

# load snRNA data
sce <- readRDS("sce_dge_input.rds")
sce
pb <- readRDS("sce_dge_pb_input.rds")
pb


###################################################################### 
#   R--> Sample-level analysis: Pseudobulk methods                   #
######################################################################
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts(
    Low_vs_Not_AD = Low-Not_AD,
    Inter_vs_Not_AD = Inter-Not_AD,
    High_vs_Not_AD = High-Not_AD,
    Inter_vs_Low = Inter-Low,
    High_vs_Inter = High-Inter,
    High_vs_Low = High-Low,
    levels = mm)

# run DS analysis
res = pbDS(pb, design = mm, contrast = contrast)
names(res)
# # run DS analysis
# #res <- pbDS(pb, verbose = FALSE)
saveRDS(res,file="sample_level_pb_dge_pair_out.rds")



# Results filtering    eg.   Low-Not_AD                                  
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.5)
  dplyr::arrange(u, p_adj.loc)
})
saveRDS(tbl_fil,file="sample_level_pb_dge_pair_Low-Not_AD_filter_out.rds")

