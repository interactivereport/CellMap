#' Estimate the cell type proportions of mixture bulk RNA
#' 
#' This function estimates cell type proportions of mixture bulk RNA samples based on one sc/sn RNAseq expression matrix.
#' 
#' @param strBulk The full path to the query mixture bulk expression file. 
#'   Expression matrix separated by tabs with rows are genes, columns are samples. 
#'   First row includes the sample names or ids, while first column consists of gene symbols or gene ensembl id.
#' @param strSC The paths to the expression matrix of a sc/sn RNAseq datasets (.rds).
#'   The expression matrix with rows are genes (official gene symbol is required as first column);
#'   and columns are cells with cell type and sample information encoded in to the column names (cellType|sampleID|…).
#'   The sample information will be used by MuSiC.
#' @param strPrefix The prefix with path of the result files. 
#'   There are two files produced: a pdf file contains all cell type decomposition figures; 
#'   a tab separated table file including composition and p-values.
#' @param method 'CellMap' or 'MuSiC' to be choosen for decomposition. Default is CellMap.
#' @param geneNameReady A boolean to indicate if the gene names in the query mixture bulk expression matrix
#'   is official symbol already.
#'   The \code{FALSE} option also works with the official symbol is used in the expression matrix.
#'   Default is \code{FALSE}, which enable to find official symbol by an R package called \code{biomaRt}.
#' @param ensemblPath The path to a folder where ensembl gene definition is/will be saved.
#'   The ensembl gene definition file will be saved if it never run before.
#'   Default is \emph{Data/} in the current working directory.
#' @param ensemblV The version of the ensembl to be used for the input query bulk expression. Default is 97.
#' @param expCutoff The minimal total measure of a cell type to be considered. Default is 1M
#' @param cellMapDEG One from ‘edgeR’, ‘DESeq2’, ‘voom’ or ‘Top’ can be chosen.
#'   This indicates the method for identifying the cell type signature genes.
#'   Since it is from one dataset, and no batch correction, default is DESeq2.
#' @param cellMapGeneN A numeric indicates the maximin number of signature genes for each cell type in an iteration.
#'   Default is 100.
#' @param core The number of computation nodes could be used. Default is 2.
#' 
#' @examples
#' strBulk <- system.file("extdata","bulk.txt",package="cellmap")
#' strSC <- system.file("extdata","GSE103723.rds",package="cellmap")
#' cellmapOne(strBulk,strSC,strPrefix="~/bulk_decomp",core=16)
#' cellmapOne(strBulk,strSC,strPrefix="~/bulk_decomp",method="MuSiC",core=16)
#' 
#' 
#' @export
cellmapOne <- function(strBulk,
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
    strSC <- paste("profile/",strSC,".rds",sep="")
    if(!file.exists(strSC)) stop(paste0(oriSC," cannot be found!"))
  }
  scID <- gsub("\\.rds","",basename(strSC))
  strF <- paste(strPrefix,"_",method,"_",scID,".pdf",sep="")

  ## decomposition -----
  BiocParallel::register(BiocParallel::MulticoreParam(core))
  if(method=="MuSiC"){
    SC <- oneSC(strSC,expCutoff)
    print(system.time({cellComp <- deconv.MuSiC(X,list(SC=SC))}))
  }else if(method=="CellMap"){
    SC <- oneBulk(strSC,expCutoff)
    selG <- featureSel(SC,
                       base::sapply(strsplit(colnames(SC),"\\|"),head,1),
                       method=cellMapDEG,
                       para=list(DEGlogFCcut=1,
                                 DEGqvalcut=0.05,
                                 DEGbasemeancut=16,
                                 selFeatureN=cellMapGeneN))
    selGn <- ceiling(min(sapply(selG,length))/100)*100
    selG <- base::sapply(selG,function(x)return(x[1:min(selGn,length(x))]))
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
    cellCol <- setNames(RColorBrewer::brewer.pal(length(cName),'Set1'),cName)
  }else if(length(cName)<13){
    cellCol <- setNames(RColorBrewer::brewer.pal(length(cName),'Set3'),cName)
  }else{
    cellCol <- setNames(scales::hue_pal()(length(cName)),cName)
  }
  pdf(strF,width=5+round(ncol(X)/3))
  plotComposition(cellComp$composition,
                  cellCol,pV=cellComp$compoP,
                  ggplotFun=ggplot2::ggtitle(paste(method," on ",scID," feature coverage:",cellComp$coverR,"%",sep="")))
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




