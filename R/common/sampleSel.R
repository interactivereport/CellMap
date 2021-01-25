########################
## sampleSel.R
##
########################

sampleSel <- function(Mix,features,para){
  mixBulk <- Mix$mBulk
  mixR <- Mix$mixR
  rm(Mix)
  ## training ---------
  traRes <- deconv(mixBulk,features,perm=para$setN,modelForm=para$modelForm)
  compos <- traRes$rawComp
  selSets <- traRes$rawSets
  #saveRDS(list(preCom=compos,expCom=mixR),file=paste(para$strfix,"rds",sep="."))
  rm(traRes)
  ## preprocessing ----
  evalRes <- selSampleID <- c()
  for(bID in names(compos)){
    estComp <- compos[[bID]]$prop
    estComp <- apply(estComp,2,function(x){return(x/sum(x))})
    estComp[compos[[bID]]$pV>0.05] <- 0
    knowR <- setNames(rep(0,nrow(estComp)),rownames(estComp))
    knowR[rownames(mixR)] <- mixR[,bID]/sum(mixR[,bID])
    res <- evalEstimation(estComp,knowR)
    ## select the good ones
    ix <- head(order(res),para$topN)#,decreasing=T
    print(summary(res[ix]))
    selSampleID <- rbind(selSampleID,selSets[ix,])
    barD <- rbind(reshape2::melt(as.matrix(estComp)[,ix])[,c(1,3)],
                  data.frame(value=res[ix],Var1=rep("Eva",length(ix))))
    ## keep training results
    predictR <- apply(estComp[,ix],1,median)
    predictR <- predictR/sum(predictR)
    evalRes <- rbind(evalRes,c(predictR-knowR,Cor=median(res[ix])))
  }
  rownames(evalRes) <- names(compos)

  return(list(selS=selSampleID,evalDiff=evalRes))
}

evalEstimation <- function(est,truR){
  w <- rep(1,length(truR))
  w[truR==0] <- 2## penalize the false positive detection
  return(apply(est,2,function(x){
    w1 = w
    w1[x==0] <- 2 #penalize the false negative
    return(sqrt(mean((w1*x-truR)^2)))
    #return(1-dist(rbind(x*w,truR))/length(truR))
  }))
  #return(cor(est,truR))
}
  
  