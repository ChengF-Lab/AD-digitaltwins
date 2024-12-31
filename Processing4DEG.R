
######################################################
## Step1: h5ad to csv (python)
######################################################

## Load pacakge and data
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
# import anndata2ri
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import sparse
import h5py
sc.settings.verbosity = 3  #verbosity:errors (0),warnings(1),info (2),hints (3)
sc.logging.print_header()


adata = sc.read_h5ad("adata_4_R_obj.h5ad")

np.max(adata.X)
t=adata.X.toarray()
pd_t = pd.DataFrame(data=t, index=adata.obs['sample_id'], columns=adata.var['gene_ids'])
pd_t.to_csv('expression.csv.gz',compression='gzip')

## metadata
meta_dat = pd.DataFrame()
meta_dat['sample_id']= adata.obs['sample_id']
meta_dat['ADNC']= adata.obs['ADNC']
meta_dat['celltype2']= adata.obs['celltype2']
meta_dat['celltype']= adata.obs['celltype']
meta_dat['donor_id']= adata.obs['Donor ID']
meta_dat.to_csv('meta.csv',index=False)

######################################################
## Step2: R--> Transfer csv to SCE   
######################################################
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

raw_dat <- fread("expression.csv.gz")
raw_meta <- read.csv("meta.csv")
class(raw_dat)
class(raw_meta)
dim(raw_dat)
dim(raw_meta)
head(raw_meta)
raw_dat[1:5,1:5]
row.names(raw_dat)=raw_dat$sample_id
names = row.names(raw_dat)
raw_dat = as.matrix(raw_dat,rownames=1)

### create sce
counts <- as.matrix(t(raw_dat))
colnames(counts) = names
counts[1:5,1:5]
dim(counts)
#[1]   36601 1378211
class(counts)
sce <- SingleCellExperiment(counts)
sce
#sce$sample_id=as.factor(paste0(raw_meta$ABC_score_class,sep="_", raw_meta$donor_id))
sce$donor_id = as.factor(raw_meta$donor_id)
sce$ADNC = raw_meta$ADNC
sce$celltype2 = raw_meta$celltype2
sce$celltype = raw_meta$celltype
names(assays(sce))="counts"
sce
#rawsce = sce
saveRDS(sce,"snRNA_sce.rds")

######################################################
## Step3: R--> filter SCE data  
######################################################

# load snRNA data
seRNA <- readRDS("snRNA_sce.rds")
seRNA

# View the first few rows of the column data to find the cell type column
head(colData(seRNA))
levels(factor(colData(seRNA)$ADNC))


colData(seRNA)$ADNC <- gsub("Intermediate", "Inter", colData(seRNA)$ADNC)
# Verify the change
table(colData(seRNA)$ADNC)


### # remove undetected genes
dim(seRNA) # [1]   35418 1191153
seRNA <- seRNA[rowSums(counts(seRNA) > 0) > 0, ]
dim(seRNA) 

qc <- perCellQCMetrics(seRNA)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
seRNA <- seRNA[, !ol]
dim(seRNA) 

# remove lowly expressed genes
seRNA <- seRNA[rowSums(counts(seRNA) > 1) >= 10, ]
dim(seRNA) 

saveRDS(seRNA,file= "snRNA_filter_sce.rds")

######################################################
## Step4: R--> Prepare Input for muscat (DEG analysis) 
######################################################
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

sce <- readRDS("snRNA_filter_sce.rds")

# ### prepare dge input
sce$id <- paste0(sce$ADNC,sce$donor_id)
(sce <- prepSCE(sce,
    kid = "celltype2", # subpopulation assignments
    gid = "ADNC",  # group IDs (ctrl/stim)
    sid = "id",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))
nk;ns;ng
names(kids) <- kids; names(sids) <- sids

# # overview the data
# # nb. of cells per cluster-sample
head(t(table(sce$cluster_id, sce$sample_id)))


saveRDS(sce,"sce_dge_input.rds")

######################################################################
#   Step6: R--> Aggregation of sc to pseudobulk data                 #
######################################################################
# ## preapre pb input
pb <- aggregateData(sce,
    assay = "counts", fun = "sum",
    by = c("cluster_id", "sample_id"))

# # one sheet per subpopulation
assayNames(pb)
t(head(assay(pb)))

# ## save data
saveRDS(pb,"sce_dge_pb_input.rds")


# ### Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))

ggsave("pbMDS_plot1.pdf")

# # use very distinctive shaping of groups & change cluster colors
pb_mds <- pb_mds +
  scale_shape_manual(values = c(17, 4, 15, 8)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# change point size & alpha
pb_mds$layers[[1]]$aes_params$size <- 5
pb_mds$layers[[1]]$aes_params$alpha <- 0.6
pb_mds
ggsave("pbMDS_plot2.pdf")



