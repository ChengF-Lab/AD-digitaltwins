#####################################
# Functions for GO term analysis
#####################################


suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
  library(stringr)
})


##########################################################################################
# 2. Load function
##########################################################################################


# Functions for creating 'low-overlapping aggregates' of cells

computeKNN <- function(data=NULL, query=NULL, k=50, includeSelf=FALSE, ...){
  # Compute KNN for query points (usually a reduced dims matrix)
  # This returns a matrix of indices mapping query to neighbors in data
  # If query has n cells (rows) and k = 50, will be a n x 50 matrix
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


getLowOverlapAggregates <- function(proj, target.agg=500, k=100, overlapCutoff=0.8, dimReduc="IterativeLSI", seed=1){
  # Generate low-overlapping aggregates of cells
  ##############################################
  # proj = ArchR project
  # target.agg = number of target aggregates (before filtering)
  # k = number of cells per aggreagate
  # overlapCutoff = Maximum allowable overlap between aggregates
  set.seed(seed)

  # Get reduced dims:
  rD <- getReducedDims(proj, reducedDims=dimReduc)

  # Subsample
  idx <- sample(seq_len(nrow(rD)), target.agg, replace = !nrow(rD) >= target.agg)

  # Get KNN Matrix:
  knnObj <- computeKNN(data=rD, query=rD[idx,], k=k)

  # Check whether aggregates pass overlap cutoff
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]

  # Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  # Name aggregates and return as a df of cell ids x aggs
  names(knnObj) <- paste0("agg", seq_len(length(knnObj)))
  knnDF <- data.frame(knnObj)[,c(3,2)]
  colnames(knnDF) <- c("cell_name", "group")
  knnDF$cell_name <- as.character(knnDF$cell_name)
  knnDF
}

#####################################
# Matrix Correlation tools
#####################################

# Rcpp script for performing fast row-wise pearson correlations:

sourceCpp(code='
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Borrowed from https://github.com/GreenleafLab/ArchR/blob/master/src/Correlation.cpp
// who in turn adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
// [[Rcpp::export]]
Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {

  if(X.ncol() != Y.ncol()){
    stop("Columns of Matrix X and Y must be equal length!");
  }

  if(max(idxX) > X.nrow()){
    stop("Idx X greater than nrow of Matrix X");
  }

  if(max(idxY) > Y.nrow()){
    stop("Idx Y greater than nrow of Matrix Y");
  }

  // Transpose Matrices
  X = transpose(X);
  Y = transpose(Y);

  const int nx = X.ncol();
  const int ny = Y.ncol();

  // Centering the matrices
  for (int j = 0; j < nx; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }

  for (int j = 0; j < ny; ++j) {
    Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
  }

  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss_X(nx);
  for (int i = 0; i < nx; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
  }

  Rcpp::NumericVector inv_sqrt_ss_Y(ny);
  for (int i = 0; i < ny; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
  }

  //Calculate Correlations
  const int n = idxX.size();
  Rcpp::NumericVector cor(n);
  for(int k = 0; k < n; k++){
    cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X( idxX[k] - 1) * inv_sqrt_ss_Y( idxY[k] - 1);
  }

  return(cor);

}'
)


cor2Matrices <- function(mat1, mat2, subset1=NULL, subset2=NULL){
    # Calculate row correlations and associated statistics of two distinct matrices using Rcpp
    ###########################
    # mat1: first matrix with rows to be correlated
    # mat2: second matrix with rows to be correlated
    # subset1: vector of rownames or indices to use for calculating correlations
    # subset2: vector or rownames or indices to use for calculating correlations

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat1) == 0)){
      stop("Error: matrix 1 contains rows of all 0's!")
    }
    if(any(Matrix::rowSums(mat2) == 0)){
      stop("Error: matrix 2 contains rows of all 0's!")
    }
    if(is.null(rownames(mat1))){
        rownames(mat1) <- 1:nrow(mat1)
    }
    mat1 <- as.matrix(mat1)
    if(is.null(rownames(mat2))){
        rownames(mat2) <- 1:nrow(mat2)
    }
    mat2 <- as.matrix(mat2)

    if(is.null(subset1)){
      subset1 <- rownames(mat1)
    }
    if(is.null(subset2)){
      subset2 <- rownames(mat2)
    }
    mat1 <- mat1[subset1,]
    mat2 <- mat2[subset2,]

    # Get indices to correlate
    message("Determining indices to correlate...")
    idx <- expand.grid(rownames(mat1), rownames(mat2))
    idx1 <- match(idx[,1], rownames(mat1))
    idx2 <- match(idx[,2], rownames(mat2))

    df <- data.frame(
        "x" = rownames(mat1)[idx1],
        "y" = rownames(mat2)[idx2]
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(idx1, idx2, mat1, mat2)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat1)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat1) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("x", "y", "Correlation", "FDR")]
    return(df)
}


getFreqs <- function(x){
  # Return a named vector of frequencies of x
  tab <- table(x) %>% as.data.frame.table()
  frqs <- tab$Freq
  names(frqs) <- tab[,1]
  frqs[order(frqs, decreasing=TRUE)]
}

theme_BOR <- function(base_size=14, base_family="Helvetica", border = TRUE) {
  library(grid)
  library(ggthemes)
  # Should plots have a bounding border?
  if(border){
    panel.border <- element_rect(fill = NA, color = "black", size = 0.7)
    axis.line <- element_blank()
  }else{
    panel.border <- element_blank()
    axis.line <- element_line(color = "black", size = 0.5)
  }

  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = panel.border,
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = axis.line,
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))

}


cmaps_BOR <- list(
  # Many of these adapted from ArchR ColorPalettes.R by Jeff Granja or colors.R from BuenColors
  # https://github.com/GreenleafLab/ArchR/blob/master/R/ColorPalettes.R
  # https://github.com/caleblareau/BuenColors/blob/master/R/colors.R

  ## Sequential colormaps:
  solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC',
                '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'), #buencolors
  sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
              "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  horizon = c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5',
              '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),
  horizonExtra =c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5",
                  "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"),
  blueYellow = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                  '#A2E700','#FFFF00','#FFD200','#FFA500'), #buencolors
  wolfgang_basic = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4",
                    "#1D91C0", "#225EA8", "#253494", "#081D58"), #buencolors
  wolfgang_extra = c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF",
                    "#2A9EC1", "#206AAD", "#243996", "#081D58"), #buencolors
  whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b'),
  whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                '#3690c0','#0570b0','#045a8d','#023858'),
  whiteViolet = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1',
                  '#DD3497', '#AE017E', '#7A0177', '#49006A'),
  comet = c("#E6E7E8","#3A97FF","#8816A7","black"),

  flame_flame = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6',
    '#FE629D', '#FF9B64', '#FFD52C', '#FFFF5F'), # buencolors

  flame_short = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6',
    '#FE629D', '#FF9B64', '#FFD52C'), # Stop short of yellow (better for tracks, etc.)

  #7-colors
  greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                '#0868ac','#084081'),

  #6-colors
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),

  #5-colors
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
  greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
  fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A"),
  purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"),
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
  zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"), #wesanderson
  darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), #wesanderson
  rushmore = c("#E1BD6D", "#EABE94", "#0B775E","#35274A" , "#F2300F"), #wesanderson
  FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"), #wesanderson
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"), #wesanderson
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"), #wesanderson
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),

  # Divergent sequential:
  coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"),
  brewer_yes = c("#053061", "#2971B1", "#6AACD0","#C1DDEB", "#F7F7F7",
                  "#FACDB5", "#E58267", "#BB2933", "#67001F"), #buencolors
  brewer_celsius = c("#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF",
                      "#FDD384", "#F88D51", "#DE3F2E", "#A50026"), #buencolors
  flame_blind = c("#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF",
                  "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"), #buencolors
  solar_flare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FFFFFF',
                  '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'), #buencolors
  brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7',
                '#FACDB5', '#E58267', '#BB2933', '#67001F'), #buencolors

  ## Qualitative colormaps:

  # see: https://carto.com/carto-colors/
  cartoPrism = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
                  '#008695', '#CF1C90', '#F97B72', '#4B4B8F'),
  cartoSafe = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99',
                 '#999933', '#882255', '#661100', '#6699CC'),
  cartoBold = c('#7F3C8D' ,'#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
                 '#008695', '#CF1C90', '#f97b72', '#4b4b8f'),
  cartoAntique = c('#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C',
                    '#9C9C5E', '#A06177', '#8C785D', '#467378'),
  cartoPastel = c('#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F', '#9EB9F3', '#FE88B1',
                   '#C9DB74', '#8BE0A4', '#B497E7', '#D3B484'),
  cartoVivid = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B',
                  '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E'),
  # 15 color
  circus = c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288",
              "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),
  iron_man = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',"#F88D51","#FAD510","#FFFF5F",'#88CFA4',
               '#238B45',"#02401B","#0AD7D3","#046C9A", "#A2A475", 'grey35'),
  # The following 3 were designed by Ryan Corces.
  stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
               "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416", "#E6C2DC"),
  calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
           "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"),
  kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
            "#FF7A5C", "#53377A", "#FF8E00","#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13")
)

calcTopGo <- function(
    allGenes, interestingGenes=NULL, pvals=NULL, geneSel=NULL,
    nodeSize=5, ontology="BP",
    alg="weight01", stat="fisher", topNodes=50
    ){
    # Calculate GO term enrichments using topGO on provided data
    # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
    ############################################################
    # allGenes: vector of genenames to be used in GO term search. Expects gene 'symbol'
    # interestingGenes: predefined list of 'instersting' genes. Incompatible with supplying pvalues.
    # geneSel: function for selecting 'interesting' genes. Can only really be a p-value cutoff...
    # pvals: vector of pvalues corresponding to geneList. If not provided, will assign everything to 1
    # nodeSize: will prune terms that have less than nodeSize number of genes
    # ontology: which GO ontology to use (MF, BP, CC)
    # alg: algorithm to be used for testing GO terms (topGO default is 'weight01')
    # stat: test statistic to use for significant GO terms
    # topNodes: how many GO terms to return in result table

    # Prepare geneList as expected for topGO (i.e. value vector with names of genes)
    if(!is.null(interestingGenes)){
        message(sprintf("Running GO enrichments with %s genes in universe of %s...",
            length(interestingGenes), length(allGenes)))
        geneList <- factor(as.integer(allGenes %in% interestingGenes))
        names(geneList) <- allGenes
        # Create topGOdata object
        GOdata <- suppressMessages(new(
            "topGOdata",
            ontology = ontology,
            allGenes = geneList,
            annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
            nodeSize = nodeSize
            ))
    }else{
        geneList <- pvals
        names(geneList) <- allGenes
        message(sprintf("Running GO enrichments with %s genes in universe of %s...",
            sum(geneSel(geneList)), length(allGenes)))
        GOdata <- suppressMessages(new(
            "topGOdata",
            ontology = ontology,
            allGenes = geneList,
            geneSel = geneSel,
            annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
            nodeSize = nodeSize
            ))
    }

    # Test for enrichment using Fisher's Exact Test
    GOresult <- suppressMessages(runTest(GOdata, algorithm=alg, statistic=stat))
    GenTable(GOdata, pvalue=GOresult, topNodes=topNodes, numChar=1000)
}


topGObarPlot <- function(goRes, cmap = NULL, nterms=10, border_color="black",
    barwidth=0.85, title="", enrichLimits=c(0.0, 5.5), barLimits=NULL){
    # Plot GO results in bar plot form
    goRes$log2FoldEnrichment <- log2(goRes$Significant / goRes$Expected)
    goRes$log2FoldEnrichment <- ifelse(goRes$log2FoldEnrichment > enrichLimits[2], enrichLimits[2], goRes$log2FoldEnrichment)
    goRes$threshPval <- ifelse(goRes$pvalue == "< 1e-30", 1e-30, as.numeric(goRes$pvalue))
    goRes$log10pval <- -log10(goRes$threshPval)
    if(!is.null(barLimits)){
        goRes$log10pval <- ifelse(goRes$log10pval < barLimits[2], goRes$log10pval, barLimits[2])
    }

    # Only plot the top nterms (reverse order to plot most significant at top)
    goRes <- goRes[1:nterms,]
    goRes <- goRes[nrow(goRes):1,]

    if(is.null(cmap)){
        cmap <- cmaps_BOR$comet
    }
    p <- (
        ggplot(goRes, aes(x=Term, y=log10pval, fill=log2FoldEnrichment))
        + geom_bar(stat="identity", width=barwidth, color=border_color)
        + scale_x_discrete(
            limits=goRes$Term, # Required to prevent alphabetical sorting of terms
            labels= function(x) str_wrap(x, width=40) # wrap long GO term labels
            )
        + scale_fill_gradientn(colors=cmap, limits=enrichLimits)
        + xlab("")
        + ylab("-log10 pvalue")
        + ggtitle(title)
        + theme_BOR(border=FALSE)
        + theme(panel.grid.major=element_blank(),
            panel.grid.minor= element_blank(),
            plot.margin = unit(c(0.25,1,0.25,1), "cm"),
            #aspect.ratio = 6/nterms, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1))
        + coord_flip()
    )
    if(!is.null(barLimits)){
        p <- p + scale_y_continuous(limits=barLimits, expand=c(0,0))
    }else{
        p <- p + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
    }
    p
}

