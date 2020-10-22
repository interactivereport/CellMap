#####################
## cellMapTraining.R
##
## main function to create cell type profiles
#####################
loadCellMapTraining <- function(){
  require(BiocParallel)
  require(RColorBrewer)
  require(colorspace)
  require(sva)
  require(edgeR)
  require(DESeq2)
  require(nnls)

  #if(!require(ggpubr)||!require(BiocParallel)||!require(optparse)) stop("missing packages!")
  source("common/createProfile.R")
  source("common/mergeProfile.R")
  source("common/evalProfile.R")
  source("common/updateProfile.R")
  source("common/plotCorProfile.R")
}
cellMapTraining <- function(strData,
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
                            rmseCutoff = 0.1, # remove the profile combination by RMSE performance
                            tailR = 0.8, # the ratio of profile sets for each sample to be remove in real Bulk evaluation
                            maxRMrate=0.75, ## maxinum ratio of total profile sets can be removed from real bulk evaluation
                            
                            maxIteration = 10, #max iteration for the full training
                            core=2
                            ){
  suppressMessages(suppressWarnings(loadCellMapTraining()))
  version <- "v1.0.5"
  para <- sapply(as.list(match.call())[-1],eval)
  register(MulticoreParam(para$core))
  para$trainD <- para$strData
  
  if(is.null(cellTypeMap)) message("The cell type name mapping is not provided, ALL cell type names from data file will be used directly!")
  para$trainMap <- para$cellMap <- cellTypeMap
  iter <- 0
  bestProfile <- allProfile <- NULL
  while(iter<para$maxIteration){
    iter <- iter + 1

    ## obtain the initial profile
    D <- initProfile(para)
    colnames(D$expr) <- paste(colnames(D$expr),"|iter",iter,sep="")
    cType <- unique(sapply(strsplit(colnames(D$expr),"\\|"),head,1))
    if(is.null(para$cellMap)) para$cellMap <- setNames(cType,cType)
    if(is.null(para$cellCol)){
      if(length(cType)<10){
        para$cellCol <- setNames(RColorBrewer::brewer.pal(length(cType),'Set1'),cType)
      }else{
        para$cellCol <- setNames(colorspace::rainbow_hcl(length(cType)),cType)
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
    ## save the profile
    strRDS <- paste0(strPrefix,"_CellMapProfile.",version,".rds")
    pdf(gsub("rds$","pdf",strRDS))
    a <- plotCorProfile(bestProfile,modelForm=para$modelForm)
    a <- dev.off()
    saveRDS(bestProfile,file=strRDS)
    ## 
    cat("Finished iteration",iter,"\n\n\n\n\n")
    if(length(para$cellMap)<1) break
  }
    
}



