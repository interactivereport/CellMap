#####################
## cellMapTraining.R
##
## main function to create cell type profiles
#####################
source("common/createProfile.R")
source("common/mergeProfile.R")
source("common/evalProfile.R")
source("common/updateProfile.R")
source("common/plotCorProfile.R")
source("common/evalCellMap.R")

source("common/pseudoBulk.R")
source("common/featureSel.R")
source("common/sampleSel.R")
source("common/batchRM.R")

source("common/deconv.R")
source("common/plotProfileDist.R")

source("common/condensGene.R")
source("common/plotEval.R")

cellmapTraining <- function(strData,
                            strPrefix,
                            cellTypeMap=NULL,
                            cellCol=NULL,
                            modelForm="log2",
                            DEGmethod="edgeR", # options: edgeR, DESeq2, Top
                            batchMethod="Full", # options: Full, None, Partial and Separate, combat_seq
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
                            
                            maxIteration = 30, #max iteration for the full training
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



