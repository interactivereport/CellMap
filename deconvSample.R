#!/usr/bin/env Rscript

## initialization of env ------------
loadDeconvSample <- function(){
  if(!require(ggpubr,quietly=T, warn.conflicts=F)) install.packages("ggpubr",repos="https://cloud.r-project.org/")
  if(!require(BiocParallel,quietly=T, warn.conflicts=F)) install.packages("BiocParallel",repos="https://cloud.r-project.org/")
  if(!require(optparse,quietly=T, warn.conflicts=F)) install.packages("optparse",repos="https://cloud.r-project.org/")
  if(!require(ggpubr)||!require(BiocParallel)||!require(optparse)) stop("missing packages!")
}
suppressMessages(loadDeconvSample())
rm(list=ls())
graphics.off()
closeAllConnections()

bDebug <- F
wd <- getwd()
initial.options <- commandArgs(trailingOnly = FALSE)
exePath <- "./"#dirname(gsub("--file=","",grep("file=",initial.options,value=T)))#"
setwd(exePath)
source("common/deconv.R")
source("common/condensGene.R")
source("common/plotComposition.R")
setwd(wd)
## get input parameters ----------------
opls <- list()

opls[[length(opls)+1]] <- make_option(c('--mixBulk','-b'), action="store", type="character", dest="strBulk", default=NULL, help="The path to the file (txt:'Tab' separated) which contains mix bulk expression (row:genes; column: samples). [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--feature','-f'), action="store", type="character", dest="strFeature", default=paste(exePath,"/Data/NNLS.features.rds",sep=""), help="The full path to the feature file (rds/txt: 'tab' separated) for the deconvolution method. [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", default="", help="prefix of the result file name [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--delCT','-d'), action="store", type="character", dest="delCT", default=NULL, help="Path to write output files to [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=4, help="Path to write output files to [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "deconvSample.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)
if(is.null(args$strBulk)){
  args$strBulk <- "Data/IHC/Data/IHC.txt"#"/camhpc/ngs/projects/GSE73721hg38/rsem_expected_count.txt"#"biogen/TST11345/genes.estcount_table.txt"
  args$strFeature <- "Training/SCDC.rds"#"Training/NNLS.linear.D_F.iter5.rds"
  args$outputPath <- "Data/IHC/"#dirname(args$strBulk)
  args$core <- 8#32
}
if(is.null(args$strBulk)) {
  message("Incomplete arguments list!!!!!!!")
  print_help(opts)
  q()
}
cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"),sep="")
cat("\n")
cat("\n")

if(nchar(args$prefix)>1)
  args$prefix <- paste(args$prefix,"_",sep="")
if(!is.null(args$delCT))
  args$delCT <- gsub(" ","",unlist(strsplit(args$delCT,",")))
modelForm <- sapply(strsplit(basename(args$strFeature),"\\."),head,1)
register(MulticoreParam(args$core))
## process features -----------------
if(grepl("rds$",args$strFeature)){
  features <- readRDS(args$strFeature)
  #features$para$decovMethod
  modelForm <- features$para$modelForm
}else{
  features <- read.table(args$strFeature,sep="\t",header=T,as.is=T,row.names=1)
  cName <- unique(sapply(strsplit(colnames(features),"\\|"),head,1))
  cellCol <- setNames(scales::hue_pal()(length(cName)),cName)
  features <- list(fM=features,para=list(cellCol=cellCol))
}

## deconvolution --------
X <- condensGene(read.table(args$strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names=1),86)
strF <- paste(args$outputPath,"/",args$prefix,features$para$decovMethod,sep="")
if(features$para$decovMethod=="NNLS") strF <- paste(strF,"_",modelForm,"_",features$para$rmBatch,sep="")#,"_",features$para$selFeature
strF <- paste(strF,".deconv.pdf",sep="")
pdf(gsub("profile.pdf","pdf",strF),width=12)
runTime <- system.time({cellComp <- deconv(X,features,features$para$decovMethod,rmCellType=args$delCT,modelForm=modelForm)})
if(!is.null(cellComp$profileDist)) for(one in cellComp$profileDist) print(one)
a <- dev.off()
if(bDebug) saveRDS(cellComp,file=gsub("pdf","rds",strF))

#conn <- file(gsub("pdf$","command",strF),"w")
#cat(paste0(names(args),": ",args,"\n"),sep="",file=conn)
#cat("\n\n")
#cat(paste0(names(runTime),": ",runTime,"\n"),sep="",file=conn)
#close(conn)

pdf(strF,width=5+round(ncol(X)/3))
plotComposition(cellComp$composition,
                features$para$cellCol,pV=cellComp$compoP,
                ggplotFun=ggtitle(paste(features$para$decovMethod,": feature coverage:",cellComp$coverR,"%",sep="")))
tmp <- dev.off()
#write.csv(cellComp$composition,file=gsub("pdf$","csv",strF))
conn <- file(gsub("pdf$","txt",strF),"w")
cat(paste('#deconvSample.R',paste(commandArgs(trailingOnly = TRUE),collapse=" ")),"\n",sep="",file=conn)
cat(paste(features$para$decovMethod,": feature coverage:",cellComp$coverR,"%",sep=""),"\n",file=conn)
write.table(apply(cellComp$composition,2,function(x){return(x/sum(x))}),file=conn,sep="\t",col.names = NA,quote=F)
if(!is.null(cellComp$compoP)){
  cat("pValue:\n",file=conn)
  write.table(cellComp$compoP,file=conn,sep="\t",col.names = NA,quote=F)
}
if(!is.null(cellComp$overallP)){
  cat("overall p-value:\n",file=conn)
  write.table(matrix(cellComp$overallP,nrow=1,dimnames=list("overallP",names(cellComp$overallP))),
              file=conn,sep="\t",col.names = NA,quote=F)
}
cat("\n\nmissing decompotion genes:\n",paste(cellComp$missingF,collapse="\n"),file=conn,sep="")
close(conn)
message("Decomposition on all samples is finished successfully!")

