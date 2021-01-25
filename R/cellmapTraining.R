#' Create a cellmap profile for desired cell types from sc/sn RNAseq data.
#' 
#' This function creates a cellmap profile include specified cell types from a set of sc/sn RNAseq data.
#'   The TPM of full length or counts of 3’end sc/sn RNAseq data is recommended.
#' 
#' @param strData A vector of paths to the expression matrix of all sc/sn RNAseq datasets (.rds).
#'   Each expression matrix with rows are genes (official gene symbol is required as first column);
#'   and columns are cells with cell type and data set information encoded in to the column names (cellType|dataset|…).
#' @param strPrefix A string indicates the prefix with path of the result files. 
#'   There are three files produced: two PDF files contains figures of the profile quality as well as
#'   performance on pseudo mixture and input bulk if provided;
#'   and an RDS file contains the profile which can be provided to cellmap function.
#' @param cellTypeMap A named vector indicates the cell types of the profile which the cellmap needed to train for.
#'   The names of the vector are the cell type names defined in the column names from expression matrix of RDS files,
#'   while the values are the cell type will be used in the final profile. 
#'   For instance, many Exhibitory and Inhibitory cells are both defined in the data,
#'   while the neuron is one of the interested cell types. 
#'   Thus, we can create a vector c(Exhibitory=Neuron, Inhibitory=Neuron) to 
#'   let the cellmap know all Exhibitory and Inhibitory cells are now called Neuron. 
#'   If NULL, all cell types defined in the data matrix will be used as original name. 
#'   Please note that: Neuron and Neurons will be considered different cell types. Default is NULL.
#' @param cellCol A named vector indicate the color of the cell types. 
#'   The names of the vector are the cell type names while the value are R-color, “#FFFFFF” is preferred. 
#'   If NULL, colors will be assigned to each of cell types. Default is NULL.
#' @param modelForm ‘linear’ or ‘log2’ for the profile model.
#'   log2 is preferred since some genes might have dominant expression values which will bias towards those genes. 
#'   Default is log2.
#' @param DEGmethod One from ‘edgeR’, ‘DESeq2’, ‘voom’ or ‘Top’ can be chosen.
#'   This indicates the method for identifying the cell type signature genes. Default is edgeR.
#' @param batchMethod One from ‘Full’, ‘Partial’, ‘Separate’ or ‘None’ can be chosen.
#'   This indicates the method for batch removal. (Please check the publication for details.)
#'   In short, if the cell types of interested are mostly overlapped among datasets, ‘Full’ is preferred, 
#'   while ‘Separate’ is for minimal overlap. Default is ‘Full’.
#' @param sampleN A numeric indicates the number of pseudo pure samples to be generated for each cell type from each dataset.
#'   Default is 5.
#' @param seqDepth A numeric indicates the total measurements (counts) for a pseudo pure sample.
#'   Default is 2M.
#' @param normDepth A numeric indicates the sequence depth to normalize the pseudo pure samples, such as CPM.
#'   Default is 1M.
#' @param geneCutoffCPM A numeric indicates the minimal normalized expression of a gene to be considered.
#'   Default is 4.
#' @param geneCutoffDetectionRatio A numeric indicates the minimal ratio of data sets where a gene expressed for a cell type.
#'   Default is 0.8.
#' @param selFeatureN A numeric indicates the maximin number of signature genes for each cell type in an iteration.
#'   Default is 100.
#' @param DEGlogFCcut A numeric indicates the minimal log fold change for a gene to be considered signature. Default is 1.
#' @param DEGqvalcut A numeric indicates the maximin FDR for a gene to be considered signature. Default is 0.05.
#' @param DEGbasemeancut A numeric indicates the minimal average expression for ‘DESeq2’, ‘voom’ or ‘Top’ methods. Default is 16.
#' @param mixN A numeric indicates the number of pseudo mixture to be generated for each dataset during pseudo training process.
#'   Default is 10.
#' @param setN A numeric indicates the number of random combinations of pseudo pure sample during pseudo training process
#'   for each iteration. Default is 1000.
#' @param topN A numeric indicates the number of the top performing (based on RMSE) combinations of pseudo pure samples
#'   for each pseudo mixture are kept during the pseudo training. Default is 50.
#' @param strBulk The path to the expression matrix with known cell type compositions.
#'   Expression matrix is tab separated with genes in rows, samples in columns. Default is ‘’.
#' @param strRate The path to the cell type compositions of the above expression matrix.
#'   The composition matrix is tab separated with cell types in rows, samples in columns. Default is ‘’.
#' @param geneNameReady A boolean to indicate if the gene names in the above bulk expression
#'   matrix is official symbol already. The \code{FALSE} option also works with the official symbol is
#'   used in the expression matrix. Default is \code{FALSE}, which enable to find official symbol by 
#'   an R package called \code{biomaRt}. Default is \code{FALSE}
#' @param ensemblPath The path to a folder where ensembl gene definition is/will be saved.
#'   The ensembl gene definition file will be saved if it never run before.
#'   Default is \emph{Data/} in the current working directory.
#' @param ensemblV The version of the ensembl to be used for the input query bulk expression. Default is 97.
#' @param rmseCutoff A numeric indicates the maximin RMSE for a bulk sample.
#'   All cell types from those bulk samples whose RMSE is larger than this value will be included in 
#'   the next iteration training. Default is 0.1.
#' @param trailR The maximin ratio of pseudo pure combinations to be removed for each bulk sample during 
#'   bulk training process. Default is 0.8.
#' @param maxRMrate The maximin ratio of total pseudo pure combinations to be removed during bulk 
#'   training process. Default is 0.75.
#' @param maxIteration The maximin iterations.
#'   If the RMSE for all bulk samples are less than indicated rmseCutoff,
#'   the iteration will stopped. Default is 10.
#' @param core The number of computation nodes could be used. Default is 2.
#' @examples
#' strData <- c(system.file("extdata","GSE103723.rds",package="cellmap"),
#'   system.file("extdata","GSE104276.rds",package="cellmap"))
#' cellTypeMap <- c(Macrophage=Microglia,Astrocytes=Astrocytes,Oligodendrocytes=Oligodendrocytes,
#'   GABAergic=Neuron,Glutamatergic=Neuron,Excitatory=Neuron,Inhibitory=Neuron)
#' cellmapTraining(strData,strPrefix="~/cellmap_profile",cellTypeMap,core=16)
#' 
#' @export
cellmapTraining <- function(strData,
                            strPrefix,
                            cellTypeMap=NULL,
                            cellCol=NULL,
                            modelForm="log2",
                            DEGmethod="edgeR", # options: edgeR, DESeq2, Top
                            batchMethod="Full", # options: Full, None, Partial and Separate
                            sampleN = 5, # the number of pure bulk extracted from each cell type each data set
                            seqDepth = 2e6, # the number of sequencing depth requirement for pseudo bulk
                            normDepth = 1e6, # the number of sequencing depth for normalization (CPM)
                            geneCutoffCPM = 4, # the cutoff CPM of a gene
                            geneCutoffDetectionRatio=0.8,# the minimal ratio of samples within at least one cell type to express the gene
                            selFeatureN = 100, # the number of features to be selected for each cell type
                            DEGlogFCcut=1, # the absolute logFC cutoff to be considered DEGs (cell type genes)
                            DEGqvalcut=0.05, # the absolute logFC cutoff to be considered DEGs (cell type genes)
                            DEGbasemeancut=16, # the absolute logFC cutoff to be considered DEGs (cell type genes)
                            
                            
                            mixN = 10, # the number of mixture samples generating per dataset in training
                            setN = 1000, # the number of pure profile compinations from pseudo mix bulk in training
                            topN = 50, # the number of top profile sets to be select for each pseudo mix bulk in training
                            
                            strBulk = "", # the expression of real bulk samples with known cell type composition (separated by "\t"), first column is the gene names, first row is the sample names
                            strRate = "", # the composition of bulk samples (separated by "\t"), first column is the cell type names, first row is the sample names matching with bulk expression
                            geneNameReady=F,
                            ensemblPath="Data/",
                            ensemblV=97,
                            rmseCutoff = 0.1, # remove the profile combination by RMSE performance
                            tailR = 0.8, # the ratio of profile sets for each sample to be remove in real Bulk evaluation
                            maxRMrate=0.75, ## maxinum ratio of total profile sets can be removed from real bulk evaluation
                            
                            maxIteration = 10, #max iteration for the full training
                            core=2
                            ){
  eval(strData)
  eval(strPrefix)
  version <- "v1.0.5"
  
  para <- sapply(formals()[c(-1,-2)],eval)
  newP <- sapply(as.list(match.call())[-1],eval)
  for(i in names(newP)) para[[i]] <- newP[[i]]
  BiocParallel::register(BiocParallel::MulticoreParam(para$core))
  para$trainD <- para$strData
  para$version <- version
  
  if(is.null(cellTypeMap)) message("The cell type name mapping is not provided, ALL cell type names from data file will be used directly!")
  para$trainMap <- para$cellMap <- cellTypeMap
  iter <- 0
  bestProfile <- allProfile <- NULL
  strRDS <- paste0(strPrefix,"_CellMapProfile",".rds")#,version
  while(iter<para$maxIteration){
    iter <- iter + 1

    ## obtain the initial profile
    D <- initProfile(para)
    colnames(D$expr) <- paste(colnames(D$expr),"|iter",iter,sep="")
    cType <- unique(base::sapply(strsplit(colnames(D$expr),"\\|"),head,1))
    if(is.null(para$trainMap)) para$trainMap <- para$cellMap <- setNames(cType,cType)
    if(is.null(para$cellCol)){
      if(length(cType)<10){
        para$cellCol <- setNames(RColorBrewer::brewer.pal(length(cType),'Set1'),cType)
      }else{
        para$cellCol <- setNames(scales::hue_pal()(length(cType)),cType)
      }      
    }
    Profile <- c(D,selProfile(D,para))
    
    ## merging with previous profile
    allProfile <- mergeProfile(bestProfile,Profile)
    rm(Profile)
    gc()
    ## evaluate  the real bulk with expected ratio
    if(file.exists(para$strBulk) && file.exists(para$strRate)){
      evalV <- evalProfile(allProfile,para)
      ## updating the profile sets and training --------------
      if(is.null(bestProfile) || bestProfile$score>evalV$score){
        bestProfile <- updateProfile(allProfile,evalV,para)
      }
      bestProfile$para$score <- c(bestProfile$para$score,evalV$score)
    }else{
      bestProfile <- allProfile
      bestProfile$para=para
    }
    
    rm(allProfile)
    para <- bestProfile$para
    ## 
    cat("Finished iteration",iter,"\n\n\n\n\n")
    if(length(para$cellMap)<1) break
  }
  ## save the profile
  pdf(gsub("rds$","pdf",strRDS))
  a <- plotCorProfile(bestProfile,modelForm=para$modelForm)
  a <- dev.off()
  saveRDS(bestProfile,file=strRDS)
  pdf(gsub("rds$","evaluation.pdf",strRDS))
  evalCellMap(bestProfile)
  a <- dev.off()
  message("\n\nThe CellMap training was completed successfully!")
}



