###############################
## createProfile.NNLS.R
##
###############################
source("common/pseudoBulk.R")
source("common/featureSel.R")
source("common/sampleSel.R")
source("common/batchRM.R")

createProfile <- function(para,postfix){
  D <- initProfile(para)
  colnames(D$expr) <- paste(colnames(D$expr),"|iter",postfix,sep="")
  return(list(expr=D$expr,selG=D$selG,sets=selProfile(D,para)))
}

initProfile <- function(para){
  ## initial extraction of pure pseudo bulk ----
  cat("create Profile: initial extraction of pure pseudo bulk\n")
  resPure <- pseudoPure(setNames(para$tech[,1],rownames(para$tech)),para$sampleN,para$trainMap,
                        para$seqDepth,para$geneCutoffDetectionRatio)
  pseudoDepth <- apply(resPure$express,2,sum)
  message("\n\nfiltering samples with less than ",para$normDepth/10," reads")
  print(pseudoDepth[pseudoDepth<para$normDepth/10])
  resPure$express <- resPure$express[,pseudoDepth>para$normDepth/10]

  resCPM <- pseudoLogCPM(resPure$express,para$normDepth,
                         cutoffCPM=para$geneCutoffCPM,
                         cutoffRatio=para$geneCutoffSampleRatio)
  ratio <- apply(resPure$express[rownames(resCPM$logCPM),],2,sum)/apply(resPure$express,2,sum)
  message("\n\nfiltering samples whose CPM remove more than ",(1-para$geneCutoffDetectionRatio)*100,"%")
  print(ratio[ratio<para$geneCutoffDetectionRatio])
  selID <- names(ratio)[ratio>para$geneCutoffDetectionRatio]
  resPure$express <- resPure$express[rownames(resCPM$logCPM),selID]
  resCPM$logCPM <- resCPM$logCPM[,selID]
  resCPM$pheno <- resCPM$pheno[selID,]
  
  #para$geneList <<- unique(c(para$geneList,rownames(resCPM$logCPM)))
  ## remove the batch effects ---------
  if(para$rmBatch=="Full"){
    cat("create Profile: removing the batch effects\n")
    rmBatchD <- batchRM(resCPM$logCPM,"cData",resCPM$pheno,para$addBatchInfo)
  }else if(para$rmBatch=="None"){
    rmBatchD <- resCPM$logCPM
  }else if(para$rmBatch=="Partial"){
    cat("create Profile: removing the batch effects for some cell types\n")
    rmBatchD <- batchRM(resCPM$logCPM,"cData",resCPM$pheno,para$addBatchInfo)
    
    dcPair <- table(resCPM$pheno[!duplicated(resCPM$pheno),"cData"])
    selX <- resCPM$pheno$cData %in% names(dcPair)[dcPair==1]
    rmBatchD[,selX] <- resCPM$logCPM[,selX]
  }else if(para$rmBatch=="Separate"){
    cat("create Profile: removing the batch effects for each cell type separately\n")
    rmBatchD <- resCPM$logCPM
    for(i in unique(resCPM$pheno$cType)){
      cat("\n",i,"\n")
      ix <- resCPM$pheno$cType==i
      rmBatchD[,ix] <- batchRM(resCPM$logCPM[,ix],"cData",resCPM$pheno[ix,])
    }
  }else{
    stop("UNKNOWN model type:",para$rmBatch,"!\n")
  }
  ## Obtain the feature selection----
  cat("create Profile: Obtaining the feature selection\n")
  if(para$selFeature=="Top"||para$selFeature=="edgeR"){
    selG <- featureSel(rmBatchD,
                       resCPM$pheno$cType,
                       method=para$selFeature,
                       resCPM$pheno$cData,
                       para$selFeatureN)
  }else{
    selG <- featureSel(resPure$express,
                       resCPM$pheno$cType,
                       method=para$selFeature,
                       resCPM$pheno$cData,
                       para$selFeatureN)
  }

  return(list(expr=rmBatchD,selG=selG))
}

selProfile <- function(D,para){
  ## initial extraction of mix pseudo bulk ----
  message("create Profile: initial extraction of mix pseudo bulk")
  mRes <- pseudoMix(para$trainD,para$sampleN,para$cellMap,para$seqDepth,para$reCondens)
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
  return(selSets[!duplicated(selSets),])
}

