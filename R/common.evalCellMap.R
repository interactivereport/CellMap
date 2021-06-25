#####################
## evalCellMap.R
##
## plot the evaluation figures: correlation by cell type, RMSE for pseudo bulk and real bulk if available
###################

evalCellMap <- function(profile){
  para <- profile$para
  
  ## performance on pseudo bulk
  mRes <- pseudoMix(para$trainD,10,para$trainMap,para$seqDepth)
  expR <- estR <- c()
  for(i in names(mRes)){
    comp <- deconv(mRes[[i]]$mBulk,profile,modelForm=para$modelForm)
    r <- as.matrix(apply(comp$composition,2,function(x)return(x/sum(x))))
    r[comp$compoP>0.05] <- 0
    estR <- cbind(estR,r)
    
    r <- merge(data.frame(row.names=names(para$cellCol)),mRes[[i]]$mixR,by='row.names',all=T)
    rownames(r) <- r[,1]
    r <- r[,-1]
    r[is.na(r)] <- 0
    r <- as.matrix(apply(r,2,function(x)return(x/sum(x))))
    expR <- cbind(expR,r)
  }
  #save(estR,expR,file="tmp/eval.rdata")
  plotEval(estR,expR,para$cellCol,T)
  
  ## performance on real bulk
  if(nchar(para$strBulk)>3){
    bulk <- condensGene(read.table(para$strBulk,header=T,sep="\t",check.names=F,as.is=T,row.names = 1),
                        para$geneNameReady,ensemblV=para$ensemblV,ensemblPath=para$ensemblPath)
    comp <- deconv(bulk,profile,modelForm=para$modelForm)
    estR <- as.matrix(apply(comp$composition,2,function(x)return(x/sum(x))))
    estR[comp$compoP>0.05] <- 0
    
    expR <- as.matrix(read.table(para$strRate,header=T,sep="\t",check.names=F,as.is=T,row.names=1))[rownames(estR),]
    plotEval(estR,expR,para$cellCol)
  }
  
}
