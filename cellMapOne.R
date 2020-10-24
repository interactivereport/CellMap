###########################
## cellMapOne.R
##
## function to decomposite bulk based on ONLY one sn/sc data set
############################

loadCellMapOne <- function(){
  require(BiocParallel)
  require(RColorBrewer)
  require(colorspace)
  require(optparse)
  require(nnls)
  require(edgeR)
  require(reshape2)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)#BiocManager::install("DESeq2")
  require(biomaRt)#BiocManager::install("biomaRt")
    
  source("common/MuSiC.R")
  source("common/deconv.R")
  source("common/pseudoSC.R")
  source("common/featureSel.R")
  source("common/plotComposition.R")
  source("common/condensGene.R")
}

cellMapOne <- function(strBulk,
                       strSC,
                       strPrefix=substr(strBulk,1,nchar(strBulk)-4),
                       method="CellMap",
                       geneNameReady=F,
                       ensemblPath="Data/",
                       ensemblV=97,
                       expCutoff=1e6,
                       cellMapDEG="DESeq2",
                       cellMapGeneN=100,
                       core=2){
  eval(strBulk)
  eval(strSC)
  suppressMessages(suppressWarnings(loadCellMapOne()))
  ## input para
  a <- sapply(as.list(match.call())[-1],eval)
  strCommand <- paste(paste(names(a),sapply(a,paste,collapse=" "),sep=":"),collapse="; ")
  message(strCommand)
  ## obtain the bulk data for decomposing -----
  X <- condensGene(read.table(strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names=1),
                   geneNameReady,ensemblV=ensemblV,ensemblPath=ensemblPath)
  ## obtain the sc/sn Data ------
  oriSC <- strSC
  if(!file.exists(strSC)){
    strSC <- paste("Data/",strSC,".rds",sep="")
    if(!file.exists(strSC)) stop(paste0(oriSC," cannot be found!"))
  }
  scID <- gsub("\\.rds","",basename(strSC))
  strF <- paste(strPrefix,"_",method,"_",scID,".pdf",sep="")

  ## decomposition -----
  register(MulticoreParam(core))
  if(method=="MuSiC"){
    SC <- oneSC(strSC,expCutoff)
    print(system.time({cellComp <- deconv.MuSiC(X,list(SC=SC))}))
  }else if(method=="CellMap"){
    SC <- oneBulk(strSC,expCutoff)
    selG <- featureSel(SC,
                       sapply(strsplit(colnames(SC),"\\|"),head,1),
                       method=cellMapDEG,
                       para=list(DEGlogFCcut=1,
                                 DEGqvalcut=0.05,
                                 DEGbasemeancut=16,
                                 selFeatureN=cellMapGeneN))
    selGn <- ceiling(min(sapply(selG,length))/100)*100
    selG <- sapply(selG,function(x)return(x[1:min(selGn,length(x))]))
    write.csv(as.data.frame(lapply(selG,'length<-',selGn)),
              file=gsub("pdf","csv",strF),
              quote=F)
    SC <- log2(1+apply(SC,2,function(x)return(1e6/sum(x)*x)))
    print(system.time({cellComp <- deconv(X,list(expr=SC,selG=selG),perm=400)}))
  }else{
    message(method," is an unknown method!")
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
                  ggplotFun=ggtitle(paste(method," on ",scID," feature coverage:",cellComp$coverR,"%",sep="")))
  tmp <- dev.off()

  conn <- file(gsub("pdf$","txt",strF),"w")
  cat("#",strCommand,"\n",file=conn)
  cat(paste(method," on ",scID," feature coverage:",cellComp$coverR,"%",sep=""),"\n",file=conn)
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
  
}




