############################
## pseudoBulk.R
## total column number smaller than 20 will be considered bulk
############################

oneSC <- function(strData,seqD=1e6){

  cat("Starting",strData,"\n")
  X <- readRDS(strData)
  #X <- X[,apply(X,2,sum)>200]
  cType <- base::sapply(strsplit(colnames(X),"\\|"),head,1)
  selType <- rep(F,ncol(X))
  for(i in unique(cType)) if(sum(X[,cType==i])>seqD) selType <- selType | cType==i
  X <- X[,selType]
  if(length(unique(cType[selType]))<2) return()

  ## remove MT genes and Ribosomal genes and miRNA----
  X <- X[!grepl("^MT|^RP|^MIR",rownames(X)),order(colnames(X))]
  print(table(cType[selType]))
  return(X)
}

oneBulk <-function(strData,seqD=1e6,sampleN=10){
  
  cat("Starting",strData,"\n")
  X <- readRDS(strData)
  cType <- base::sapply(strsplit(colnames(X),"\\|"),head,1)
  selType <- rep(F,ncol(X))
  for(i in unique(cType)) if(sum(X[,cType==i])>seqD) selType <- selType | cType==i
  X <- X[,selType]
  if(length(unique(cType[selType]))<2){
    message(paste("No cell type cluster read number is over",seqD))
    return()
  }
  
  cType <- cType[selType]
  bulk <- c()
  for(j in unique(cType)){
    iPOS <- grep(paste("^",j,"$",sep=""),cType)
    if(sum(X[,iPOS])<seqD) next
    cat("\nworking on",j,sum(cType==j),"\t")
    ix <- rep(0,length(iPOS))
    one <- c()
    if(sum(X[iPOS])*2<seqD){
      message("Too few cells!")
      next
    }
    for(k in 1:sampleN){
      ix <- sample(iPOS,length(ix),replace=T)
      while(sum(X[,ix])>10*seqD && length(ix)>0.5*length(iPOS))
        ix <- sample(iPOS,length(ix)-ceiling(length(iPOS)/20),replace=T)
      while(sum(X[,ix])<seqD && length(ix)<2*length(iPOS))
        ix <- sample(iPOS,length(ix)+ceiling(length(iPOS)/10),replace=T)
      cat(length(ix),";",sep="")
      one <- cbind(one,apply(X[,ix],1,sum))
    }
    colnames(one) <- paste(j,1:sampleN,sep="|")
    bulk <- cbind(bulk,one)
  }
  rm(X)
  message("")
  ## remove MT genes and Ribosomal genes and miRNA----
  bulk <- bulk[!grepl("^MT|^RP|^MIR",rownames(bulk)),]
  ## ------
  return(bulk)
}
