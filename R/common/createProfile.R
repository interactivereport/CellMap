###############################
## createProfile.NNLS.R
##
###############################

createProfile <- function(para,postfix){
  D <- initProfile(para)
  colnames(D$expr) <- paste(colnames(D$expr),"|iter",postfix,sep="")
  if(is.null(para$cellMap)){
    cType <- unique(base::sapply(strsplit(colnames(D$expr),"\\|"),head,1))
    para$cellMap <- setNames(cType,cType)
  }
  return(list(expr=D$expr,selG=D$selG,sets=selProfile(D,para)))
}

initProfile <- function(para){
  ## initial extraction of pure pseudo bulk ----
  cat("create Profile: initial extraction of pure pseudo bulk\n")
  resPure <- pseudoPure(para$strData,para$sampleN,para$trainMap,
                        para$seqDepth,para$geneCutoffDetectionRatio)
  pseudoDepth <- apply(resPure$express,2,sum)
  message("\n\nfiltering samples with less than ",para$normDepth/10," reads")
  print(pseudoDepth[pseudoDepth<para$normDepth/10])
  resPure$express <- resPure$express[,pseudoDepth>para$normDepth/10]

  resCPM <- pseudoLogCPM(resPure$express,para$normDepth,
                         cutoffCPM=para$geneCutoffCPM,
                         cutoffRatio=para$geneCutoffDetectionRatio)
  ratio <- apply(resPure$express[rownames(resCPM$logCPM),],2,sum)/apply(resPure$express,2,sum)
  message("\n\nfiltering samples whose CPM removal is more than ",(1-para$geneCutoffDetectionRatio)*100,"%")
  print(ratio[ratio<para$geneCutoffDetectionRatio])
  selID <- names(ratio)[ratio>para$geneCutoffDetectionRatio]
  resPure$express <- resPure$express[rownames(resCPM$logCPM),selID]
  resCPM$logCPM <- resCPM$logCPM[,selID]
  resCPM$pheno <- resCPM$pheno[selID,]
  
  ## remove the batch effects ---------
  if(para$batchMethod=="Full"){
    cat("create Profile: removing the batch effects\n")
    rmBatchD <- batchRM(resCPM$logCPM,"cData",resCPM$pheno,para$addBatchInfo)
  }else if(para$batchMethod=="None"){
    rmBatchD <- resCPM$logCPM
  }else if(para$batchMethod=="Partial"){
    cat("create Profile: removing the batch effects for some cell types\n")
    rmBatchD <- batchRM(resCPM$logCPM,"cData",resCPM$pheno,para$addBatchInfo)
    
    dcPair <- table(resCPM$pheno[!duplicated(resCPM$pheno),"cData"])
    selX <- resCPM$pheno$cData %in% names(dcPair)[dcPair==1]
    rmBatchD[,selX] <- resCPM$logCPM[,selX]
  }else if(para$batchMethod=="Separate"){
    cat("create Profile: removing the batch effects for each cell type separately\n")
    rmBatchD <- resCPM$logCPM
    for(i in unique(resCPM$pheno$cType)){
      cat("\n",i,"\n")
      ix <- resCPM$pheno$cType==i
      rmBatchD[,ix] <- batchRM(resCPM$logCPM[,ix],"cData",resCPM$pheno[ix,])
    }
  }else{
    stop("UNKNOWN model type:",para$batchMethod,"!\n")
  }
  ## Obtain the feature selection----
  cat("create Profile: Obtaining the feature selection\n")
  if(para$DEGmethod=="Top"||para$DEGmethod=="edgeR"){
    para$DEGbasemeancut <- log2(1+para$DEGbasemeancut)
    selG <- featureSel(rmBatchD,
                       resCPM$pheno$cType,
                       method=para$DEGmethod,
                       resCPM$pheno$cData,
                       para)
  }else{
    selG <- featureSel(resPure$express,
                       resCPM$pheno$cType,
                       method=para$DEGmethod,
                       resCPM$pheno$cData,
                       para)
  }

  return(list(expr=rmBatchD,selG=selG))
}

selProfile <- function(D,para){
  ## initial extraction of mix pseudo bulk ----
  message("create Profile: initial extraction of mix pseudo bulk")
  mRes <- pseudoMix(para$trainD,para$mixN,para$cellMap,para$seqDepth)
  if(length(mRes)<1 || sum(sapply(mRes,function(x)return(ncol(x$mBulk))))<para$sampleN) stop("selProfile error in mixture extraction")
  ## selecting prifles -------
  message("create Profile: selecting profiles")
  trainEval <- selSets <- c()
  for(i in names(mRes)){
    message("create Profile: training ",i)
    trainRes <- sampleSel(mRes[[i]],D,para)
    selSets <- rbind(selSets,trainRes$selS)
    trainEval <- rbind(trainEval,trainRes$evalDiff)
    message("Finished Profile: training ",i," with ",nrow(trainRes$selS)," sets for ",ncol(trainRes$selS)," cell types")
    rm(trainRes)
  }
  #saveRDS(trainEval,file=paste(para$strfix,".training.rds",sep=""))
  return(list(sets=selSets[!duplicated(selSets),]))
}

