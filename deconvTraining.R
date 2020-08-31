#!/usr/bin/env Rscript
loadDeconvTraining <- function(){
  #if(!require(ggpubr,quietly=T, warn.conflicts=F)) install.packages("ggpubr",repos="https://cloud.r-project.org/")
  if(!require(BiocParallel,quietly=T, warn.conflicts=F)){
    install.packages("BiocParallel",repos="https://cloud.r-project.org/")
    if(!require(BiocParallel,quietly=T, warn.conflicts=F)) stop("BiocParallel cannot be installed!")
  }
  if(!require(RColorBrewer,quietly=T, warn.conflicts=F)){
    install.packages("RColorBrewer",repos="https://cloud.r-project.org/")
    if(!require(RColorBrewer,quietly=T, warn.conflicts=F)) stop("RColorBrewer cannot be installed!")
  }
  if(!require(colorspace,quietly=T, warn.conflicts=F)){
    install.packages("colorspace",repos="https://cloud.r-project.org/")
    if(!require(colorspace,quietly=T, warn.conflicts=F)) stop("colorspace cannot be installed!")
  }
  #if(!require(ggpubr)||!require(BiocParallel)||!require(optparse)) stop("missing packages!")
  source("common/createProfile.R")
  source("common/mergeProfile.R")
  source("common/evalProfile.R")
  source("common/updateProfile.R")
  source("common/plotCorProfile.R")
}
suppressMessages(require(optparse,quietly=T, warn.conflicts=F))
## input ----
opls <- list()
opls[[length(opls)+1]] <- make_option(c('--Data','-d'), action="store", type="character", dest="data", default=NULL,
                                      help="The full path to the file listed training sc/sn RNAseq data (required). Three columns separated by tab without header.\n
                                      \tfirst column: the name of the dataset;\n\tsecond column: the full path to the data rds file;\n\t
                                      third column:the cell types defined in the datasets (first element separated by | in column names), e.g. A, B, C,
                                      and the cell type names prepered, e.g. CT1, CT2. In format of: A=CT1;B=CT2;C=CT2.")
opls[[length(opls)+1]] <- make_option(c('--DEG','-g'), action="store", type="character", dest="selFeature", default='DESeq2',
                                      help="The method (DESeq2 or edgeR) to extract the cell type genes. [default: %default].")
opls[[length(opls)+1]] <- make_option(c('--bulk','-b'), action="store", type="character", dest="strBulk",default="",
                                      help="The full path to the real bulk expression matrix with known composition (tab separated, row is the genes, column is the samples)")
opls[[length(opls)+1]] <- make_option(c('--rate','-r'), action="store", type="character", dest="strRate",default="",
                                      help="The full path to the cell type compostion matrix (tab separated, row is the cell type, column is the samples).")
opls[[length(opls)+1]] <- make_option(c('--out','-o'), action="store", type="character", dest="output",default=NULL,
                                      help="The output CellMap feature files, including the prefix.")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=8, help="Path to write output files to [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "deconvTraining.R", description = "", epilogue = "")

para <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)
if(is.null(para$data) || is.null(para$output)) {
  message("Incomplete arguments list!")
  print_help(opts)
  q()
}
message(paste(paste(names(para),para,sep=": "),collapse="\n"))

## initial the parameters -----------
suppressMessages(loadDeconvTraining())
para$tech <- read.table(para$data,sep="\t",as.is=T,row.names=1)
print(para$tech)
## get the celltypeMap
tmp <- unlist(strsplit(para$tech[,2],";"))
cellTypeID <- sapply(strsplit(tmp,"="),head,1)
cellTypeName <- setNames(sapply(strsplit(tmp,"="),tail,1)[!duplicated(cellTypeID)],cellTypeID[!duplicated(cellTypeID)])
cType <- unique(cellTypeName)

if(length(cType)<10){
  para$cellCol <- setNames(RColorBrewer::brewer.pal(length(cType),'Set1'),cType)
}else{
  para$cellCol <- setNames(colorspace::rainbow_hcl(length(cType)),cType)
}
para$trainMap <- cellTypeName
para$trainD <- setNames(para$tech[,1],rownames(para$tech))
para$decovMethod <- 'CellMap'
para$modelForm <- 'log2'
para$rmBatch <- 'Separate'
para$seqDepth <- 2e6
para$version <- 'v1.0.2'
para <- c(para,list(
  sampleN = 5, # the number of pure bulk extracted from each cell type each data set
  normDepth = 1e6, # the number of sequencing depth for normalization (CPM)
  geneCutoffCPM = 4, # the cutoff CPM of a gene
  geneCutoffSampleRatio=0.2, # the minimal ratio of samples to express the gene
  geneCutoffDetectionRatio=0.8,
  selFeatureN = 300, # the number of features to be selected for each cell type
  addBatchInfo = NULL, # besides the data sets additional information for batch removal, "cType" or NULL
  cellMap=para$trainMap,
  iteration = 1000, # the number of permutation for testing the pseudo mix bulk
  topN = 50, # the number of top profile sets to be select in pseudo mix bulk training
  tailR = 0.8, # the ratio of profile sets to be remove in real Bulk evaluation
  maxRMrate=0.75, ## maxinum ratio of profile sets can be removed
  coeffCutoff = 0.8,
  rmseCutoff = 0.1,
  maxIteration = 10, #max iteration for the training
  geneList=c(),
  debug=F,
  reCondens=F # is the gene condensed scRNA for pure sample needed to be saved
))

register(MulticoreParam(para$core))
## CellMap training -----
iter <- 0
bestProfile <- allProfile <- NULL
while(iter<para$maxIteration&&length(para$trainMap)>0){
  iter <- iter + 1
  para$strfix <- paste(para$output,para$decovMethod,para$version,sep="_")
  strPDF <- paste(para$strfix,".pdf",sep="")
  ## obtain the initial profile
  Profile <- createProfile(para,iter)
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
  pdf(gsub("pdf$","profile.pdf",strPDF))
  a <- plotCorProfile(bestProfile,modelForm=para$modelForm)
  a <- dev.off()
  saveRDS(bestProfile,file=gsub("pdf$","rds",strPDF))
  ## 
  cat("Finished iteration",iter,"\n\n\n\n\n")
  if(length(para$cellMap)<1) break
}



