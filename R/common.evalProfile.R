###############################
## evalProfile.R
##
###############################
#source("common.deconv.R")
#source("common.condensGene.R")

evalProfile <- function(fProfile,para){
  ## deconvolution and plot-----
  cat("eval Profile: deconvolution ... \n")
  bulk <- condensGene(read.table(para$strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names = 1),
                      para$geneNameReady,ensemblV=para$ensemblV,ensemblPath=para$ensemblPath)
  Res <- deconv(bulk,fProfile,modelForm=para$modelForm)
  ## evaluation by correlation for each sample and plotting---------
  cat("eval Profile: evaluation ... \n")
  mixR <- as.matrix(read.table(para$strRate,header=T,sep="\t",check.names=F,as.is=T,row.names=1))
  tR <- Res$composition
  tR[,] <- 0
  if(sum(colnames(mixR)!=colnames(tR))>0) stop("evalProfile:inconsistency between predicted and expected sample ID")
  tR[,] <- mixR[rownames(tR),]
 
  predictR <- apply(Res$composition,2,function(x){return(x/sum(x))})#plotComposition(Res$composition,cellCol,Res$compoP,ggplotFun=ggtitle(paste("Predicted",para$version)),rateReturn=T)
  predictR[Res$compoP > 0.05] <- 0

  rmse <- getRMSE(predictR,tR)
  ## decision ----------
  cat("eval Profile: deciding ... \n")
  poorCellType <- c()
  rmSetsIndex <- c()
  
  ## using rmse ------------
  poorS <- rmse
  poorS[poorS < para$rmseCutoff] <- 0
  bugetRM <- round(poorS/sum(poorS)*para$tailR*nrow(fProfile$sets))

  for(i in names(rmse)){
    if(rmse[i]>para$rmseCutoff){
      poorCellType <- c(poorCellType,rownames(tR)[tR[,i]>0])
      rawComp <- apply(Res$rawComp[[i]]$prop,2,function(x){return(x/sum(x))})
      rawComp[Res$rawComp[[i]]$pV>0.05] <- 0
      rmSetsIndex <- c(rmSetsIndex,order(apply(rawComp,2,getRMSE,tR[,i]),decreasing=T)[1:bugetRM[i]])
    }
  }
#  saveRDS(cbind(as.data.frame(t(apply(Res$composition,2,function(x){return(x/sum(x))})-tR)),Cor=coeff),
#          file=paste(para$strfix,".eval.rds",sep=""))
  message("Poor performance on cell types:",paste(poorCellType,collapse="; "))
  return(list(score=sqrt(sum(rmse^2)*sum(rmse>para$rmseCutoff)),
              poorCellType=unique(poorCellType),
              rmSetsIndex=unique(rmSetsIndex)))
  
}
getRMSE <- function(predictR,mixR){
  if(is.null(ncol(predictR))){
    return(sqrt(mean((predictR-mixR)^2)))
  }
  rmse <- apply(mixR-predictR[rownames(mixR),colnames(mixR)],2,function(x){return(sqrt(mean(x^2)))})
  print(rmse)
  return(rmse)
}
