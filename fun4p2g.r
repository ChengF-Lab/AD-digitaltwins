# Functions for working with peak to gene linkages

getP2G_GR <- function(proj, corrCutoff=NULL, varCutoffATAC=0.25, varCutoffRNA=0.25, filtNA=TRUE){
  # Function to get peaks and genes involved in peak to gene links
  ############################################################
  # proj: ArchR project that alreayd has Peak2GeneLinks
  # corrCutoff: minimum numeric peak-to-gene correlation to return
  # varCutoffATAC: minimum variance quantile of the ATAC peak accessibility when selecting links
  # varCutoffRNA: minimum variance quantile of the RNA gene expression when selecting links
  p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
  p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
  # Remove peaks with 'NA' correlation values
  if(filtNA){
    p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
  }
  if(!is.null(corrCutoff)){
    p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
  }
  # Filter by variance quantile
  p2gDF <- p2gDF[which(p2gDF$VarQATAC > varCutoffATAC & p2gDF$VarQRNA > varCutoffRNA),]
  # The genomic range contains just the peak ranges:
  p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
  mcols(p2gGR) <- p2gDF
  p2gGR
}


grLims <- function(gr){
  # Get the minimum and maximum range from a GR
  if(length(gr) == 0){
    return(NA)
  }
  starts <- start(gr)
  ends <- end(gr)
  c(min(starts, ends), max(starts, ends))
}


getP2Gregions <- function(proj, genes, p2gGR=NULL, corrCutoff=0.4, buffer_space=0.05, min_width=25000, ...) {
  # Function to get regions containing entire peak to intertesed genes,
  # i.e. a GR that contains all peak to gene links
  ###############################################################
  # p2gGR: genomic range containing all peak to gene links
  # genes: vector of genes to look up
  # buffer_space: fraction of total length to expand on each side of region

  # Get gene GR from ArchR project
  geneGR <- promoters(getGenes(proj)) # Promoters gets 2kb upstream and 200bp downstream
  geneGR <- geneGR[!is.na(geneGR$symbol)]

  # if p2gGR not provided, pull it from ArchR project
  if(is.null(p2gGR)){
    p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff, ...)
  }

  # Now for each gene, construct GR of all loops and gene TSS
  resultGR <- geneGR[match(genes, geneGR$symbol)]
  start(resultGR) <- sapply(resultGR$symbol, function(x){
      min(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })
  end(resultGR) <- sapply(resultGR$symbol, function(x){
      max(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })

  # Finally, resize by buffer space
  resultGR <- resize(resultGR, width=width(resultGR) + buffer_space*width(resultGR), fix="center")
  resultGR <- resize(resultGR, width=ifelse(width(resultGR) > min_width, width(resultGR), min_width), fix="center")
  resultGR
}


