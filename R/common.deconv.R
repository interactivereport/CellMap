#'  deconv.R
#'  @noRd
deconv <- function(bulk,feature,perm=100,batch=NULL,rmCellType=NULL,modelForm='log2'){
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
  inTraining <- F
  if(is.null(selSets)){
    cType <- sapply(strsplit(colnames(features),"\\|"),head,1)
    cName <- sort(unique(cType))
    selSets <- c()
    for(i in cName) selSets <- cbind(selSets,colnames(features)[sample(grep(i,cType),perm,replace=T)])
    selSets <- selSets[!duplicated(selSets),]
    cat("Generated ramdom sets!\n")
    inTraining <- T
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
  bulk <- apply(bulk,2,function(x){return(1e6/sum(x)*x)})
  if(logFlag){
  	bulk <- log2(1+bulk)
  }else{
  	features <- 2^features-1
  }
  features <- as.matrix(features[comF,])
  bulk <- as.matrix(bulk[comF,])

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
  ## deconvolution ------
  cName <- base::sapply(strsplit(selSets[1,],"\\|"),head,1)
  FullProp <- BiocParallel::bplapply(data.frame(bulk,check.names=F),function(x){
    one <- base::apply(selSets,1,function(fID,x){
      M <- features[,fID]
      y <- x
      if(modelForm=="zscore"){
        y <- (x-mean(x))/sd(x)
        M <- (M-mean(M))/sd(M)
      }
      fit <- nnls::nnls(M,y)
      ## initialize the return
      prop <- setNames(fit$x,cName)
      pV <- setNames(rep(0.0499,length(cName)),cName)
      fitP <- F.pvalue(y,fit$residuals,ncol(M))
      rmse <- sqrt(sum(fit$residuals^2)/length(y))
      if(inTraining){
          res <- NLSfit(M,y,prop)
          if(!is.null(res)) return(res)
      }
      return(list(prop=prop,pV=pV,fitP=fitP,rmse=rmse))
    },x)
    return(list(prop=sapply(one,function(x){return(x$prop)}),
                pV=sapply(one,function(x){return(x$pV)}),
                fitP=sapply(one,function(x){return(x$fitP)}),
                rmse=sapply(one,function(x){return(x$rmse)})))
  })#,BPPARAM=MulticoreParam(tasks=1)
  ## combin the results and prepare the return --------
  topN <- 5
  finalProp <- base::sapply(names(FullProp),function(one){
      ix <- order(FullProp[[one]]$rmse)[1:topN]
      for(i in ix){
          res <- NLSfit(features[,selSets[i,]],bulk[,one],
                        FullProp[[one]]$prop[,i])
          if(!is.null(res)){
              FullProp[[one]]$prop[,i] <- res$prop
              FullProp[[one]]$pV[,i] <- res$pV
              FullProp[[one]]$rmse[i] <- res$rmse
              FullProp[[one]]$fitP[i] <- res$fitP
          }
      }
      w <- median(FullProp[[one]]$rmse,na.rm=T) - FullProp[[one]]$rmse[ix]
      w <- w/sum(w)
      return(list(prop=(FullProp[[one]]$prop[,ix]%*%as.matrix(w))[,1],
                  pV=(10^-(-log10(FullProp[[one]]$pV[,ix])%*%as.matrix(w)))[,1],
                  fitP=10^-sum(-log10(FullProp[[one]]$fitP[ix])%*%w),
                  rmse=sum(FullProp[[one]]$rmse[ix]*w)))
  })
  return(list(composition=apply(finalProp,2,function(one){return(one$prop)}),
              compoP=apply(finalProp,2,function(one){return(one$pV)}),
              overallP=apply(finalProp,2,function(one){return(one$fitP)}),
              rmse=apply(finalProp,2,function(one){return(one$rmse)}),
              coverR=coverR,rawComp=FullProp,rawSets=selSets,missingF=missingF,
              missingByCellType=missingByCellType,
              rmGeneQC=rmGeneQC))
}

NLSfit <- function(M,y,x){
    cName <- names(x)
    Init <- setNames(x,paste("b",1:ncol(M),sep=""))
    formu <- paste("bulk~",paste(paste(names(Init),"*",cName,sep=""),collapse="+"),sep="")
    Data <- as.data.frame(cbind(M,y))
    colnames(Data) <- c(cName,"bulk")
    
    fitR <- try(nls(as.formula(formu),Data,Init,algorithm="port",lower=rep(0,length(Init))),silent=T)
    res <- NULL
    if("convergence"%in%names(fitR)){
        a <- summary(fitR)
        res <- list()
        res$prop <- setNames(a$coefficients[,"Estimate"],cName)
        res$pV <- setNames(a$coefficients[,"Pr(>|t|)"],cName)
        res$rmse <- sqrt(sum(a$residuals^2)/length(y))
        res$fitP <- F.pvalue(y,a$residuals,ncol(M))
    }
    return(res)
}

F.pvalue <- function(y,resid,df2){
    df1 <- 1
    n <- length(y)
    RSS1 <- sum((y-mean(y))^2)
    RSS2 <- sum(resid^2)
    fstat <- ((RSS1-RSS2)/(df2-df1))/(RSS2/(n-df2))
    return(pf(fstat,df1,df2,lower.tail=F))
}



