############################
## pseudoBulk.R
## total column number smaller than 20 will be considered bulk
############################
require(BiocParallel)
#source("common/condensGene.R")

pseudoPure <- function(accIDs,sN,cellMap,seqD,gCover=0.6,condensF=F,strPreGeneFile=NULL){
  bulk <- c()
  genes <- list()
  res <- bplapply(accIDs,onePure,sN,cellMap,seqD,condensF)
  ## merge pseudo bulk require genes be reported by certain percentage of data sets
  #return(res[sapply(res,length)>0])
  gNames <- table(unlist(sapply(res,rownames)))/length(res)
  gNames <- names(gNames)[gNames>gCover]
  #strF <- "Data/subNeuron.common.gene.rds"
  if(!is.null(strPreGeneFile)&&file.exists(strPreGeneFile)){
    gNames <- gNames[gNames%in%readRDS(strF)]
  }
  cat("gName in pseudoPure:",length(gNames),"\n")
  for(i in names(accIDs)){
    cat(i,"in pseudoPure\n")
    genes[[i]] <- rownames(res[[i]])
    totalReads <- apply(res[[i]],2,sum)
    res[[i]] <- res[[i]][rownames(res[[i]])%in%gNames,]
    cat("\t",apply(res[[i]],2,sum)/totalReads,"\n")
    if(length(bulk)==0){
      bulk <- res[[i]]
    }else{
      bulk <- merge(bulk,res[[i]],by="row.names",all=T,sort=F)
      rownames(bulk) <- bulk[,1]
      bulk <- bulk[,-1]
    }
  }
  bulk[is.na(bulk)] <- 0
  return(list(express=bulk,geneL=genes))
}

pseudoLogCPM <- function(X,normDep=1e6,grp=NULL,cutoffCPM=8,cutoffRatio=0.2){
  ## normalized to gene length 
  #X <- apply(X,2,function(x,gLen){return(x/gLen)},as.numeric(sapply(strsplit(rownames(X),"\\|"),tail,1)))
  ####
  CPM <- apply(X,2,function(x){return(normDep/sum(x)*x)})
  cName <- sapply(strsplit(colnames(CPM),"\\|"),head,1)
  selGene <- rep(F,nrow(CPM))
  for(i in unique(cName)){
    selGene <- selGene | apply(CPM[,cName==i],1,function(x){return(sum(x<cutoffCPM)/length(x))})<cutoffRatio
  }
  fLogCPM <- log2(1+CPM[selGene,])
  pheno <- data.frame(row.names=colnames(fLogCPM),
                      cType=cName,
                      cData=sapply(strsplit(colnames(fLogCPM),"\\|"),"[[",2))
  cat("pure logCPM:",range(fLogCPM),"\n")
  return(list(logCPM=fLogCPM,pheno=pheno))
}

#names of cellMap is the ones in the data, the value is the one used for doconvolution
pseudoMix <- function(accIDs,sN,cellMap,seqDep,condensF=F){ 
  mixR <- bulk <- list()
  res <- bplapply(accIDs,oneMix,sN,cellMap,seqDep,condensF)
  res[sapply(res,is.null)] <- NULL
  return(res)
}

onePure <-function(dID,sampleN,cellMap,seqD=2e6,recondense=F){
  strData <- dID#paste("Data/",dID,".rds",sep="")
  dID <- gsub("\\.rds","",basename(dID))
  if(!file.exists(strData)){
    cat("Missing",strData,"\n")
    return()
  }
  cat("Starting",strData,"\n")
  X <- readRDS(strData)
  if(recondense&&as.Date(file.mtime(strData),"%Y-%m-%d")<Sys.Date()){
    X <- condensGene(readRDS(strData))
    saveRDS(X,file=strData)
  }
  bulk <- c()
  cType <- sapply(strsplit(colnames(X),"\\|"),head,1)
  ix <- cType %in% names(cellMap)
  cType[ix] <- cellMap[cType[ix]]
  for(j in intersect(unique(cellMap),unique(cType))){
    cat("\tworking on",j,sum(cType==j),":")
    if(ncol(X)<30){# for bulk
      one <- as.matrix(X[,cType==j])
      colnames(one) <- sapply(strsplit(colnames(one),"\\|"),function(x){return(paste(c(j,x[-1],x[1]),collapse="|"))})
    }else{
      iPOS <- grep(j,cType)
      ix <- rep(0,ceiling(length(iPOS)))
      one <- c()
      for(k in 1:sampleN){
        ix <- sample(iPOS,length(ix),replace=T)
        while(sum(X[,ix])>10*seqD && length(ix)>0.5*length(iPOS))
          ix <- sample(iPOS,length(ix)-ceiling(length(iPOS)/20),replace=T)
        while(sum(X[,ix])<seqD && length(ix)<2*length(iPOS))
          ix <- sample(iPOS,length(ix)+ceiling(length(iPOS)/10),replace=T)
        cat(length(ix),"; ",sep="")
        one <- cbind(one,apply(X[,ix],1,sum))
      }
      colnames(one) <- paste(j,dID,1:sampleN,sep="|")
    }
    bulk <- cbind(bulk,one)
    cat("\n")
  }
  rm(X)
  ## remove MT genes and Ribosomal genes and miRNA----
  bulk <- bulk[!grepl("^MT|^RP|^MIR",rownames(bulk)),]
  cat("\tTotal range:",range(bulk),"\n")
  ## ------
  return(bulk)
}

extCellIndex <- function(nCell,cType){
  ix <- c()
  cTypeN <- table(cType)
  iR <- setNames(sample(100,length(cTypeN)),names(cTypeN))
  while(sum(iR/sum(iR)<0.05)>0) iR <- setNames(sample(100,length(cTypeN)),names(cTypeN))
  iR <- round(nCell * iR/sum(iR))
  for(i in names(iR)){
    ix <- c(ix,sample(grep(paste("^",i,"$",sep=""),cType),iR[i],replace=T))
  }

  ## select 10% of nCell for each type ----
#  for(i in unique(cType)) ix <- c(ix,sample(grep(i,cType),ceiling(0.1*nCell),replace=T))
  ## select the result of nCell from the total ----
#  ix <- c(ix,sample(length(cType),nCell-length(ix),replace=T))
  ## ----
  return(ix)
}
oneMix <- function(dID,sampleN,cellMap,seqD=2e6,recondense=F){
  strData <- dID#paste("Data/",dID,".rds",sep="")
  dID <- gsub("\\.rds","",basename(dID))
  if(!file.exists(strData)){
    cat("Missing",strData,"\n")
    return()
  }
  selType <- names(cellMap)
  cat("Starting extracting mixture from",dID,"\n")
  X <- readRDS(strData)
  if(recondense&&as.Date(file.mtime(strData),"%Y-%m-%d")<Sys.Date()){
    X <- condensGene(readRDS(strData))
    saveRDS(X,file=strData)
  }
  
  selIndex <- sapply(strsplit(colnames(X),"\\|"),head,1)%in%selType
  if(sum(selIndex)<10) return(NULL)
  X <- X[,selIndex]
  cType <- cellMap[sapply(strsplit(colnames(X),"\\|"),head,1)]
  print(table(cType))
  ## selecting ----
  selCell <- bulk <- c()
  nCell <- as.integer(ncol(X)*0.5)
  for(i in 1:sampleN){
    ix <- extCellIndex(nCell,cType)
    while(sum(X[,ix])<seqD){
      nCell <- nCell + as.integer(ncol(X)*0.1)
      ix <- extCellIndex(nCell,cType)
    }
    bulk <- cbind(bulk,apply(X[,ix],1,sum,na.rm=T))
    res <- table(cType[ix])
    selCell <- merge(selCell,
                     setNames(data.frame(res,row.names=names(res))[,-1,drop=F],LETTERS[i]),
                     all=T,by="row.names")
    rownames(selCell) <- selCell[,1]
    selCell <- selCell[,-1,drop=F]
  }
  rm(list=c("X"))
  selCell[is.na(selCell)] <- 0
  colnames(selCell) <- colnames(bulk) <- paste(dID,1:sampleN,sep="|")
  ## remove MT genes and Ribosomal genes and miRNA----
  bulk <- bulk[!grepl("^MT|^RP|^MIR",rownames(bulk)),]
  ## ------
  print(selCell)
  return(list(mBulk=bulk,mixR=selCell))
}




