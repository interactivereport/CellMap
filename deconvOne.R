#!/usr/bin/env Rscript
## deconvOne.R
rm(list=ls())
graphics.off()
closeAllConnections()
loadDeconvOne <- function(){
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
  if(!require(optparse,quietly=T, warn.conflicts=F)){
    install.packages("optparse",repos="https://cloud.r-project.org/")
    if(!require(optparse,quietly=T, warn.conflicts=F))  stop("optparse cannot be installed!")
  }
  source("common/MuSiC.R")
  source("common/deconv.R")
  source("common/pseudoSC.R")
  source("common/featureSel.R")
  source("common/plotComposition.R")
  source("common/condensGene.R")
  message("loading completed!")
}
suppressMessages(loadDeconvOne())
fDebug <- F
strCommand <- paste('#deconvOne.R',paste(commandArgs(trailingOnly = TRUE),collapse=" "))
## get input parameters ----------------
opls <- list()

opls[[length(opls)+1]] <- make_option(c('--mixBulk','-b'), action="store", type="character", dest="strBulk", default=NULL, help="The path to the file (txt:'Tab' separated) which contains mix bulk expression (row:genes; column: samples). [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--sc','-f'), action="store", type="character", dest="strSC", default='GSE103723',
                                      help="accession ID or path to the sc/sn RNA expression data. The sc/sn RNAseq accession ID,available: GSE76381, GSE103723, phs001836 and GSE67835; the sc/sn RNA expression data should be rds matrix, with genes as row, cell as columns (column names separated by '|' with the first element to be cell type) [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--cutoff','-c'), action="store", type="numeric", dest="cutoff", default=1e6, 
                                      help="Numeric: the total reads/measurements of all cells for a celltype to be considered [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--method','-m'), action="store", type="character", dest="method", default="MuSiC", 
                                      help="The methods to be used to decompose (MuSiC or CellMap) [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", default=paste(getwd(),"CellMap_",sep=""), help="prefix of the result file name [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=8, help="Numeric: the number of cores to be used [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "deconvOne.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

if(fDebug){
  args$strBulk <- "/home/zouyang/CellMap/biogen/TST11452/genes.tpm_table.txt"
  args$prefix <- "/home/zouyang/CellMap/biogen/TST11452"
  args$method="CellMap"
  args$core <- 1
}
if(is.null(args$strBulk)) {
  message("Incomplete arguments list!")
  print_help(opts)
  q()
}
## obtain the bulk data for decomposing -----
X <- condensGene(read.table(args$strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names=1),99)
## obtain the sc/sn Data ------
strSC <- paste("Data/",args$strSC,".rds",sep="")
if(file.exists(args$strSC)){
  strSC <- args$strSC
}
if(!file.exists(strSC)){
  message(args$strSC," cannot be found!")
  q()
}
scID <- gsub("\\.rds","",basename(strSC))
if(dir.exists(args$prefix)){
  strF <- paste(args$prefix,"/",args$method,"_",scID,".pdf",sep="")
}else{
  strF <- paste(args$prefix,"_",args$method,"_",scID,".pdf",sep="")
}

## decomposition -----
register(MulticoreParam(args$core))
if(args$method=="MuSiC"){
  SC <- oneSC(strSC,args$cutoff)
  print(system.time({cellComp <- deconv.MuSiC(X,list(SC=SC))}))

}else if(args$method=="CellMap"){
  SC <- oneBulk(strSC,args$cutoff)
  selG <- featureSel(SC,
                     sapply(strsplit(colnames(SC),"\\|"),head,1),
                     method='DESeq2',
                     featureN=300)
  selGn <- ceiling(min(sapply(selG,length))/100)*100
  selG <- sapply(selG,function(x)return(x[1:min(selGn,length(x))]))
  write.csv(as.data.frame(lapply(selG,'length<-',selGn)),
              file=gsub("pdf","csv",strF),
              quote=F)
  SC <- log2(1+apply(SC,2,function(x)return(1e6/sum(x)*x)))
  print(system.time({cellComp <- deconv(X,list(expr=SC,selG=selG),perm=400)}))
}else{
  message(args$method," is an unknown method!")
}
## save the results --------
cName <- sort(unique(sapply(strsplit(colnames(SC),"\\|"),head,1)))
if(length(cName)<9){
  cellCol <- setNames(brewer.pal(length(cName),'Set1'),cName)
}else if(length(cName)<13){
  cellCol <- setNames(brewer.pal(length(cName),'Set3'),cName)
}else{
  cellCol <- setNames(scales::hue_pal()(length(cName)),cName)
}
pdf(strF,width=5+round(ncol(X)/3))
plotComposition(cellComp$composition,
                cellCol,pV=cellComp$compoP,
                ggplotFun=ggtitle(paste(args$method," on ",scID,": feature coverage:",cellComp$coverR,"%",sep="")))
tmp <- dev.off()

conn <- file(gsub("pdf$","txt",strF),"w")
cat(strCommand,"\n",file=conn)
cat(paste(args$method," on ",scID,": feature coverage:",cellComp$coverR,"%",sep=""),"\n",file=conn)
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



