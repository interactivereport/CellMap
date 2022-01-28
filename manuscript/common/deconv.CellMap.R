##########################
## deconv.NNLS.R
##
##########################
loadNNLS <- function(){
  if(!require(nnls)){
    install.packages("nnls",repos="https://cloud.r-project.org/")
    if(!require(nnls)) stop("Cannot install nnls!")
  }
  if(!require(BiocParallel)){
    install.packages("BiocParallel",repos="https://cloud.r-project.org/")
    if(!require(BiocParallel)) stop("Cannot install BiocParallel!")
  }
  if(!require(edgeR)){
    install.packages("edgeR",repos="https://cloud.r-project.org/")
    if(!require(edgeR)) stop("Cannot install edgeR!")
  }
}
suppressWarnings(suppressMessages(loadNNLS()))
deconv.CellMap <- function(bulk,feature,perm=100,modelForm="log2",rmCellType=NULL,
                        debug=F){#c("Astrocytes","Endothelial","Macrophage","Neuron","Oligodendrocytes")
  #save(bulk,feature,file="tmp/decov.rdata")
  if(!is.null(rmCellType)){
    feature$sets <- feature$sets[,!sapply(strsplit(feature$sets[1,],"\\|"),head,1)%in%rmCellType]
    feature$selG <- feature$selG[!names(feature$selG)%in%rmCellType]
  }
  features <- feature$expr
  DEG <- feature$selG[order(names(feature$selG))]

  selG <- unique(unlist(DEG))
  selSets <- feature$sets
  cellCol <- feature$para$cellCol
  anchor <- feature$anchor
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
  rmGene <- grepl("^MT|^RP|^MIR",rownames(bulk))#|^LINC|^ATP
  rmGeneQC <- apply(bulk[rmGene,],2,sum)/apply(bulk,2,sum)
  bulk <- log2(1+apply(bulk[!rmGene,],2,function(x)return(x/sum(x)*1e6)))
  
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
  
  
  ## RUVg normalization, no normalization if anchor is NULL
  if(F){
    res <- batchRUVgAnchor(bulk,features,anchor,exec=F)
    bulk <- as.matrix(res$bulk[selG,,drop=F])
    features <- as.matrix(res$profile[selG,,drop=F])
  }else{
    bulk <- as.matrix(bulk[comF,])
    features <- as.matrix(features[comF,])
  }
  
  if(modelForm!="log2"){
    bulk <- 2^bulk -1
    features <- 2^features-1
  }
  ## add the average expression profile to the end ----- 
  if(!is.null(cellCol)){
    oneMedian <- sapply(names(cellCol),function(i){return(apply(features[,grepl(i,colnames(features))],1,median))})
    colnames(oneMedian) <- paste(colnames(oneMedian),"median",sep="|")
    features <- cbind(features,oneMedian)
    selSets <- rbind(selSets,colnames(oneMedian))
  }
  
  profileDist <- NULL
  #save(features,cellCol,DEG,file="Training/test.rdata")
  if(!is.null(cellCol)&&debug){
    source("common/plotProfileDist.R")
    for(i in names(DEG)) DEG[[i]] <- intersect(DEG[[i]],rownames(features))
    plotBulkFeature(bulk,features,cellCol,DEG)
    profileDist <- plotProfileDist(features,cellCol,DEG)
  }
  ## deconvolution ------
  cName <- sapply(strsplit(selSets[1,],"\\|"),head,1)
  FullProp <- bplapply(data.frame(bulk,check.names=F),function(x){
    one <- apply(selSets,1,function(fID,x){
      M <- features[,fID]
      y <- x
      fit <- nnls(M,y)
      ## initialize the return
      prop <- setNames(fit$x,cName)
      pV <- setNames(rep(0.0499,length(cName)),cName)
      fitP <- F.pvalue(y,fit$residuals,ncol(M))
      rmse <- sqrt(sum(fit$residuals^2)/length(y))
      return(list(prop=prop,pV=pV,fitP=fitP,rmse=rmse))
    },x)
    return(list(prop=sapply(one,function(x){return(x$prop)}),
                pV=sapply(one,function(x){return(x$pV)}),
                fitP=sapply(one,function(x){return(x$fitP)}),
                rmse=sapply(one,function(x){return(x$rmse)})))
  })#,BPPARAM=MulticoreParam(tasks=1)
  ## combin the results and prepair the reture --------
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
  
  selx <- unique(unlist(base::sapply(FullProp,function(one)return(order(one$rmse)[1:topN]))))
  heatData <- list(bulk=bulk,features=features[,sort(unique(as.vector(selSets[selx,])))])

  return(list(composition=apply(finalProp,2,function(one){return(one$prop)}),
              compoP=apply(finalProp,2,function(one){return(one$pV)}),
              overallP=apply(finalProp,2,function(one){return(one$fitP)}),
              rmse=apply(finalProp,2,function(one){return(one$rmse)}),
              coverR=coverR,rawComp=FullProp,rawSets=selSets,missingF=missingF,
              missingByCellType=missingByCellType,
              rmGeneQC=rmGeneQC,
              profileDist=profileDist,
              heatData=heatData))
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



