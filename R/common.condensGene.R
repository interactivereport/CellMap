####################
## condensGene.R
## mapping to Ensembl_gene_id, hgnc_id (HGNC:...), EntrezGene_id to HGNC_SYMBOL 
####################

condensGene <- function(X,gReady=F,fun=max,ensemblV=97,ensemblPath=NULL){
  X[is.na(X)] <- 0
  if(gReady) return(X)
  ## obtain the gene definition information from biomaRt
  attrib <- c("ensembl_gene_id","hgnc_id","entrezgene_id","entrezgene","hgnc_symbol","transcript_length","chromosome_name")
  if(!is.null(ensemblPath)){
    strF <- paste(ensemblPath,"/ensembl.v",ensemblV,".rds",sep="")
    if(file.exists(strF)){
      geneSym <- readRDS(strF)
    }else{
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensemblV)
      attrib <- attrib[attrib%in%biomaRt::listAttributes(ensembl)[,1]]
      geneSym <- biomaRt::getBM(attrib,mart = ensembl)#,"transcript_length"
      if(!dir.exists(dirname(strF))) dir.create(dirname(strF))
      saveRDS(geneSym,file=strF)
    }
  }else{
    ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensemblV)
    attrib <- attrib[attrib%in%biomaRt::listAttributes(ensembl)[,1]]
    geneSym <- biomaRt::getBM(attrib,mart = ensembl)#,"transcript_length"
  }

  ## obtain the gene definition in the matrix -------
  gID <- rownames(X)
  if(is.character(X[1,1])){
    gID <- X[,1]
    X <- X[,-1,drop=F]
  }
  cat("Obtained gene definition from data\n")
  ## remove _ and . --------
  for(removeS in c("_","\\.","\\|")){
    gID <- base::sapply(strsplit(gID,removeS),head,1)
  }
  ## match data gene definition to hgnc symbol -------
  gName <- rep(NA,length(gID))
  for(i in colnames(geneSym)){
    if(i=="transcript_length") next
    gName[is.na(gName)] <- geneSym[match(gID[is.na(gName)],geneSym[,i]),"hgnc_symbol"][1:sum(is.na(gName))]
  }
  cat("ensembl (v",ensemblV,") coverage: ",sum(!is.na(gName))," in ",length(gName),sep="","\n")
  ## remove the genes without symbol definition and duplicated gene symbols -----
  selG <- !is.na(gName) & nchar(gName)>1
  X <- X[selG,]
  gName <- gName[selG]
  cat("Removing gene duplicates from data ... \n")
  for(i in unique(gName[duplicated(gName)])){
    #cat("\t",i,"\n")
    ix <- gName == i
    if(ncol(X)==1) X[ix,1][1] <- fun(X[ix,])
    else X[ix,][1,] <- apply(X[ix,],2,fun)
  }
  X <- X[!duplicated(gName),,drop=F]
  rownames(X) <- gName[!duplicated(gName)]
  ## append max gene length to row names -----
  #rownames(X) <- paste(rownames(X),geneSym[match(rownames(X),geneSym[,3]),4],sep="|")
  ## ----
  return(X)
}

rmDuplicates <- function(X){
  gName <- rownames(X)
  cat("Removing gene duplicates from data ... \n")
  pos <- c()
  uGene <- unique(gName[duplicated(gName)])
  Y <- matrix(0,nrow=length(uGene),ncol=ncol(X),
              dimnames=list(uGene,colnames(X)))
  n <- 0
  for(i in uGene){
    ix <- grep(paste('^',i,"$",sep=""),gName)# gName == i
    n <- n+1
    if(n%%100==0){
      cat(n,"/",length(uGene),"\t")
      cat("\t",i,":",length(ix),"\n")
    }
    if(ncol(X)==1) X[ix,1][1] <- max(X[ix,])
    else{
      Y[i,] <- apply(X[ix,],2,max)
      pos <- c(pos,ix[1])
    }
  }
  if(length(pos)>0) X[pos,] <- Y
  rm(Y)
  gc()
  X <- X[!duplicated(gName),,drop=F]
  gc()
  #rownames(X) <- gName[!duplicated(gName)]
  return(X)
}


