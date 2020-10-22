##########################
## deconv.R
##
##########################

deconv <- function(bulk,feature,perm=100,batch=NULL,rmCellType=NULL,modelForm='log2',
                        TMM=F,FullG=F,debug=F){#c("Astrocytes","Endothelial","Macrophage","Neuron","Oligodendrocytes")
  if(!is.null(rmCellType)){
    feature$sets <- feature$sets[,!sapply(strsplit(feature$sets[1,],"\\|"),head,1)%in%rmCellType]
    feature$selG <- feature$selG[!names(feature$selG)%in%rmCellType]
  }
  features <- feature$expr
  DEG <- feature$selG[order(names(feature$selG))]
  
  selG <- unique(unlist(DEG))
  selSets <- feature$sets
  cellCol <- feature$para$cellCol
  rm(feature)
  if(is.null(selSets)){
    cType <- sapply(strsplit(colnames(features),"\\|"),head,1)
    cName <- sort(unique(cType))
    selSets <- c()
    for(i in cName) selSets <- cbind(selSets,colnames(features)[sample(grep(i,cType),perm,replace=T)])
    selSets <- selSets[!duplicated(selSets),]
    cat("Generated ramdom sets!\n")
  }
  cat(nrow(selSets),"sets!\n")
  ## remove the MT/RP/MIR genes
  rmGene <- grepl("^MT|^RP|^MIR",rownames(bulk))
  rmGeneQC <- apply(bulk[rmGene,],2,sum)/apply(bulk,2,sum)
  bulk <- bulk[!rmGene,]
  ## common genes----
  comF <- intersect(selG,rownames(bulk))
  coverR <- 100*round(length(comF)/length(selG),3)
  cat(coverR,"% features is covered.\n",sep="")
  if(length(comF)<50){
    cat("No prediction: Too little feature coverage!\n")
    return()
  }
  missingF <- selG[!selG%in%comF]
  missingByCellType <- list()
  for(i in names(DEG)) missingByCellType[[i]] <- DEG[[i]][!DEG[[i]]%in%comF]
  ## normalize feature profiles with input data ---
  logFlag <- F
  if(modelForm=="log2") logFlag <- T
  if(TMM){
    D <- merge(2^features-1,bulk,by="row.names",all=T,sort=F)
    rownames(D) <- D[,1]
    D <- D[,-1]
    D[is.na(D)] <- 0
    D <- cpm(calcNormFactors(DGEList(D)),log=logFlag)
    features <- as.matrix(D[comF,colnames(features)])
    bulk <- as.matrix(D[comF,colnames(bulk)])
    rm(D)
  }else{
    bulk <- apply(bulk,2,function(x){return(1e6/sum(x)*x)})
    if(logFlag){
      bulk <- log2(1+bulk)
    }else{
      features <- 2^features-1
    }
    features <- as.matrix(features[comF,])
    bulk <- as.matrix(bulk[comF,])
  }
  ## add the average expression profile to the end ----- 
  if(!is.null(cellCol)){
    oneMedian <- sapply(names(cellCol),function(i){return(apply(features[,grepl(i,colnames(features))],1,median))})
    colnames(oneMedian) <- paste(colnames(oneMedian),"median",sep="|")
    features <- cbind(features,oneMedian)
    selSets <- rbind(selSets,colnames(oneMedian))
  }
  
  ## combat normalization if batch infor is provided ---------
  if(!is.null(batch)){
    X <- cbind(features,bulk)
    X <- batchRM(X,"cData",data.frame(row.names=colnames(X),cData=batch))
    features <- X[,colnames(features)]
    bulk <- X[,colnames(bulk)]
    rm(X)
  }
  profileDist <- NULL
  #save(features,cellCol,DEG,file="Training/test.rdata")
  if(!is.null(cellCol)&&debug){
    source("common/plotProfileDist.R")
    for(i in names(DEG)) DEG[[i]] <- intersect(DEG[[i]],rownames(features))
    plotBulkFeature(bulk,features,cellCol,DEG)
    #profileDist <- plotProfileDist(features,cellCol,DEG)
  }
  ## deconvolution ------
  cName <- sapply(strsplit(selSets[1,],"\\|"),head,1)
  FullProp <- bplapply(data.frame(bulk,check.names=F),function(x){
    one <- apply(selSets,1,function(fID,x){
      M <- features[,fID]
      y <- x
      if(modelForm=="zscore"){
        y <- (x-mean(x))/sd(x)
        M <- (M-mean(M))/sd(M)
      }
      fit <- nnls(M,y)
      ## initialize the return
      prop <- setNames(fit$x,cName)
      pV <- setNames(rep(0.0499,length(cName)),cName)
      fitP <- 1
      yhat <- M %*% as.matrix(fit$x)
      rmse <- sqrt(sum((y-yhat)^2)/length(y))
      ##
      Init <- setNames(fit$x,paste("b",1:ncol(M),sep=""))
      formu <- paste("bulk~",paste(paste(names(Init),"*",cName,sep=""),collapse="+"),sep="")
      Data <- as.data.frame(cbind(M,y))
      colnames(Data) <- c(cName,"bulk")
      fitR <- try(nls(as.formula(formu),Data,Init,algorithm="port",lower=rep(0,length(Init))),silent=T)
      if("convergence"%in%names(fitR)){
        a <- summary(fitR)
        prop[] <- a$coefficients[,"Estimate"]
        pV[] <- a$coefficients[,"Pr(>|t|)"]
        yhat <- predict(fitR)
        fitP <- cor.test(yhat,y)$p.value
        rmse <- sqrt(sum((y-yhat)^2)/length(y))
      }
      return(list(prop=prop,pV=pV,fitP=fitP,rmse=rmse))
    },x)
    return(list(prop=sapply(one,function(x){return(x$prop)}),
                pV=sapply(one,function(x){return(x$pV)}),
                fitP=sapply(one,function(x){return(x$fitP)}),
                rmse=sapply(one,function(x){return(x$rmse)})))
  })#,BPPARAM=MulticoreParam(tasks=1)
  ## combin the results and prepair the reture --------
  finalProp <- base::sapply(FullProp,function(one){
    #cat(one$rmse)
    ix <- order(one$rmse)[1:5]
    w <- median(one$rmse,na.rm=T) - one$rmse[ix]
    w <- w/sum(w)
    return(list(prop=(one$prop[,ix]%*%as.matrix(w))[,1],
                pV=(10^-(-log10(one$pV[,ix])%*%as.matrix(w)))[,1],
                fitP=10^-sum(-log10(one$fitP[ix])%*%w),
                rmse=sum(one$rmse[ix]*w)))
  })
  return(list(composition=apply(finalProp,2,function(one){return(one$prop)}),
              compoP=apply(finalProp,2,function(one){return(one$pV)}),
              overallP=apply(finalProp,2,function(one){return(one$fitP)}),
              rmse=apply(finalProp,2,function(one){return(one$rmse)}),
              coverR=coverR,rawComp=FullProp,rawSets=selSets,missingF=missingF,
              missingByCellType=missingByCellType,
              rmGeneQC=rmGeneQC,
              profileDist=profileDist))
}




