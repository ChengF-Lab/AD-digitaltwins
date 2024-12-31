## The code and supporting data for paper "Single-cell digital twins identify new targets and repurposable drugs in Alzheimer’s disease"


### Code
The code contains the following data analysis steps:
- snRNA-seq processing: `snRNA_QC.py`
- Differential expressed genes analysis: `Processing4DEG.R`, `DEG.R`
- snATAC-seq processing: `snATAC_QC_Clustering.R`, `snRNA_snATAC_integration.R`
- snATAC-seq peak calling: `peakCalling.R`
- Motif annotation, Peak2GeneLinkage, Positive TF_regulators: `Peak_annotation.R`
- Target gene of PosTFs: `fun4p2g.r`, `fun4TFtargerGene.r`, and `PosTF_targets.r`

### Data Description 

**Supporting data include 8 tables:**
- S. Table 1: Metadata information of 84 donors of snRNA-seq and snATAC-seq.
- S. Table 2: All significant DEGs for all cell type across AD progression.
- S. Table 3: Enrichment score of motif TFs for cell type significant marker peaks. 
- S. Table 4: Positive TF-regulators for all cell type across AD progression. 
- S. Table 5: AD-associated genes identified by overlapping AD significant SNPs with cCREs.
- S. Table 6: AD-associated lead SNPs.
- S. Table 7: Potential targets of positive TF-regulators.
- S. Table 8: Drug enrichment results by network proximity and GSEA methods.

