###########
## batchRM.R
##
## CPM[features,samples]: log2 normalized features
## batch: the column name in pheno indicates the batch
## pheno: a data frame whose row number is the same as column number of CPM, columns are interested phenotype information
#############
if(!require(sva,quietly=T, warn.conflicts=F)){
  install.packages("sva",repos="https://cloud.r-project.org/")
  if(!require(sva,quietly=T, warn.conflicts=F)) stop("sva cannot be installed!")
}

batchRM <- function(CPM,batch,pheno,additional=NULL,method="combat",core=8){
  if(sum(colnames(pheno)%in%batch)!=1){
    print(batch)
    print(colnames(pheno))
    stop("Batch is not specific in pheno!")
  }
  
  if(method=="combat"){
    if(is.null(additional)){
      combatD <- batchComBat(CPM,as.character(pheno[,batch]),core=core)
    }else{
      combatD <- batchComBat(CPM,pheno[,batch],pheno[,additional,drop=F],core=core)
    }
  }else{
    stop("Unknown normalization method:",method)
  }
  return(combatD)
}

batchComBat <- function(eData,batch,pheno=NULL,core=8){
  cat("Starting ComBat....\n")
  print(table(batch))
  if(length(table(batch))<2) return(eData)
  selX <- apply(eData,1,function(x){return(sum(x>1))})>=0.8 * min(table(batch))
  if(!is.null(pheno)){
    combatD <- ComBat(dat=eData[selX,],batch=batch,mod=model.matrix(as.formula(paste("~",paste(colnames(pheno),collapse = "+"))),data=pheno))
  }else{
    combatD <- ComBat(dat=eData[selX,],batch=batch,par.prior=F,BPPARAM=MulticoreParam(core))#
  }
  combatD[is.na(combatD)] <- 0
  eData[selX,] <- combatD
  cat("Finishing ComBat:",range(combatD),"\n")
  eData[eData<0] <- 0
  return(eData)
}



