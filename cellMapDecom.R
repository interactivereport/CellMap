######################
## cellMapDecom.R
##
## decompose the bulk based on input cell type profiles
#######################
loadCellMapDecom <- function(){
  require(ggpubr)
  require(biomaRt)#BiocManager::install("biomaRt")
  require(BiocParallel)
  require(reshape2)
  require(ggplot2)
  require(nnls)
  
  source("common/deconv.R")
  source("common/condensGene.R")
  source("common/plotComposition.R")
  
}

cellMapDecom <- function(strBulk,
                         strProfile,
                         strPrefix=substr(strBulk,1,nchar(strBulk)-4),
                         delCT=NULL,
                         geneNameReady=F,
                         ensemblPath="Data/",
                         ensemblV=97,
                         bReturn=F,
                         core=2){
  eval(strBulk)
  eval(strProfile)
  suppressMessages(suppressWarnings(loadCellMapDecom()))

  oriSC <- strProfile
  if(!file.exists(strProfile)){
    strProfile <- paste("profile/",strProfile,".rds",sep="")
    if(!file.exists(strProfile)) stop(paste0(oriSC," cannot be found!"))
  }
  profile <- readRDS(strProfile)
  X <- condensGene(read.table(strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names=1),geneNameReady,ensemblV=ensemblV,ensemblPath=ensemblPath)
  
  ## deconvolution --------
  strF <- paste0(strPrefix,"_CellMap_",profile$para$version,".pdf")
  print(system.time({cellComp <- deconv(X,profile,rmCellType=delCT,modelForm=profile$para$modelForm)}))

  #saveRDS(cellComp,file=gsub("pdf","rds",strF))
  if(bReturn) return(cellComp)
  
  pdf(strF,width=5+round(ncol(X)/3))
  plotComposition(cellComp$composition,
                  profile$para$cellCol,pV=cellComp$compoP,
                  ggplotFun=ggtitle(paste0("CellMap(",profile$para$version,") feature coverage:",cellComp$coverR,"%")))
  tmp <- dev.off()

  conn <- file(gsub("pdf$","txt",strF),"w")
  cat(paste('#Bulk:',strBulk,"Profile:",strProfile,"delCT:",paste(delCT,collapse=";"),"\n"),file=conn)
  cat(paste0("CellMap (",profile$para$version,") feature coverage: ",cellComp$coverR,"%"),"\n",file=conn)
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









