###########
## batchRM.R
##
## X[features,samples]: log2 normalized features for ComBat, counts for ComBat_seq
## batch: the column name in pheno indicates the batch
## pheno: a data frame whose row number is the same as column number of X, columns are interested phenotype information
#############

batchRM <- function(X,batch,pheno,additional=NULL,method="combat",core=8){
  if(sum(colnames(pheno)%in%batch)!=1){
    print(batch)
    print(colnames(pheno))
    stop("Batch is not specific in pheno!")
  }
  
  if(method=="combat"){
    if(is.null(additional)){
      combatD <- batchComBat(X,as.character(pheno[,batch]),core=core)
    }else{
      combatD <- batchComBat(X,pheno[,batch],pheno[,additional,drop=F],core=core)
    }
  }else if(method=='combat_seq'){
      combatD <- batchComBat_seq(X,pheno[,batch])
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
    combatD <- sva::ComBat(dat=eData[selX,],batch=batch,mod=model.matrix(as.formula(paste("~",paste(colnames(pheno),collapse = "+"))),data=pheno))
  }else{
    combatD <- sva::ComBat(dat=eData[selX,],batch=batch,par.prior=F,BPPARAM=BiocParallel::MulticoreParam(core))#
  }
  combatD[is.na(combatD)] <- 0
  eData[selX,] <- combatD
  cat("Finishing ComBat:",range(combatD),"\n")
  eData[eData<0] <- 0
  return(eData)
}

batchComBat_seq <- function(eCounts,batch){
    message("Starting ComBat_seq ...")
    print(table(batch))
    if(length(table(batch))<2) return(eCounts)
    ## genes should be expressed in 80% of samples of any batch
    selX <- apply(eCounts,1,function(x){return(sum(x>1))})>=0.8 * min(table(batch))
    combatD <- sva::ComBat_seq(counts=as.matrix(eCounts[selX,]),batch=batch)
    combatD[combatD<0] <- 0
    cat("Finishing ComBat:",range(combatD),"and demision:",dim(combatD),"\n")
    return(combatD)
}

