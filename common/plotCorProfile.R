###########################
## plotCorProfile.R
##
##################################
if(!require(reshape2,quietly=T, warn.conflicts=F)){
  install.packages("reshape2",repos="https://cloud.r-project.org/")
  if(!require(reshape2,quietly=T, warn.conflicts=F)) stop("reshape2 cannot be installed!")
}
if(!require(ggplot2,quietly=T, warn.conflicts=F)){
  install.packages("ggplot2",repos="https://cloud.r-project.org/")
  if(!require(ggplot2,quietly=T, warn.conflicts=F)) stop("ggplot2 cannot be installed!")
}

source("common/plotProfileDist.R")
plotCorProfile <- function(Profiles,selN=2000,corMethod="pearson",main="",modelForm='log2'){#spearman
  allCor <- c()
  if(is.null(Profiles$sets)){
    cType <- sapply(strsplit(colnames(Profiles$expr),"\\|"),head,1)
    cName <- sort(unique(cType))
    selSets <- c()
    for(i in cName) selSets <- cbind(selSets,colnames(Profiles$expr)[sample(grep(i,cType),selN/2,replace=T)])
    Profiles$sets <- selSets[!duplicated(selSets),]
    cat("Generated",nrow(selSets),"ramdom sets!\n")
    rm(selSets)
  }
  Profiles$expr <- as.matrix(Profiles$expr[,order(colnames(Profiles$expr))])
  Profiles$selG <- Profiles$selG[order(names(Profiles$selG))]
  if(modelForm=="linear") Profiles$expr <- 2^Profiles$expr-1
  Profiles$sets <- Profiles$sets[!duplicated(Profiles$sets),]
  sel <- melt(diag(ncol(Profiles$sets)))$value!=1
  ix <- sample(nrow(Profiles$sets),min(selN,nrow(Profiles$sets)))
  allG <- unique(unlist(Profiles$selG))
  for(i in ix){
    #if(which(ix==i)%%50==0) cat(which(ix==i),"in",selN,"\n")
    one <- sort(Profiles$sets[i,])
    res <- cor(Profiles$expr[allG,one],Profiles$expr[allG,one],method=corMethod)
    DIST <- as.matrix(dist(t(Profiles$expr[allG,one])))
    dimnames(res) <- sapply(dimnames(res),function(x){return(list(sapply(strsplit(x,"\\|"),head,1)))})
    allCor <- rbind(allCor,cbind(melt(res)[sel,],eucDist=melt(DIST)[sel,"value"]))
  }
  Prof <- Profiles$expr[allG,unique(as.vector(Profiles$sets[ix,]))]
  
  cellCol <- Profiles$para$cellCol
  if(is.null(cellCol)) return()
  
  cName <- sapply(strsplit(colnames(Prof),"\\|"),head,1)
  DS <- c()
  for(i in sort(unique(cName))){
    #cat(i,":",sum(cName==i),"\n")
    res <- cor(Prof[,cName==i],Prof[,cName==i],method=corMethod)
    DIST <- as.matrix(dist(t(Prof[,cName==i])))
    dimnames(res) <- sapply(dimnames(res),function(x){return(list(sapply(strsplit(x,"\\|"),head,1)))})
    Data <- rbind(allCor[allCor[,1]==i,],cbind(melt(res)[melt(lower.tri(res))$value,],eucDist=melt(DIST)[melt(lower.tri(DIST))$value,"value"]))
    print(ggplot(Data,aes(Var2,value,fill=Var2))+geom_boxplot(notch=T)+
            scale_fill_manual(values=cellCol)+
            scale_y_continuous(trans='log10')+
            ggtitle(paste(main,i,sep=":"))+ylab("Correlation Coefficients")+
            theme_classic()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    print(ggplot(Data,aes(Var2,eucDist,fill=Var2))+geom_boxplot(notch=T)+
            scale_fill_manual(values=cellCol)+
            scale_y_continuous(trans='log10')+
            ggtitle(paste(main,i,sep=":"))+ylab("Euclidean Distance")+
            theme_classic()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    DS <- setNames(c(DS,popDis(Data)),c(names(DS),i))
  }
  print(ggplot(cbind(as.data.frame(DS),method="NNLS",cType=names(DS)),
               aes(x=method,y=DS))+
          geom_boxplot(outlier.shape = NA)+
          geom_point(aes(color=cType),position=position_jitterdodge(jitter.width=0.2),size=2)+
          scale_color_manual(values=cellCol)+
          ylab("DS")+
          theme_classic())
  for(one in plotProfileDist(Profiles$expr[allG,],cellCol,selG=Profiles$selG)) print(one)
  return(DS)
}

popDis <- function(Data){
  cType <- as.character(Data$Var1[1])
  one <- summary(Data[Data$Var2==cType,"value"])
  mu <- one[3]
  sig <- max(one[3]-one[2],one[5]-one[3])
  D <- c()
  for(i in levels(Data$Var2)){
    if(i==cType) next
    one <- summary(Data[Data$Var2==i,"value"])
    muT <- one[3]
    sigT <- max(one[3]-one[2],one[5]-one[3])
    D <- setNames(c(D,(mu-muT)^2/(sig^2+sigT^2)),c(names(D),i))
  }
  return(min(D))
}
