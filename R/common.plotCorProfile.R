###########################
## plotCorProfile.R
##
##################################

plotCorProfile <- function(Profiles,selN=2000,corMethod="pearson",main="",modelForm='log2'){#spearman
  allCor <- c()
  if(is.null(Profiles$sets)){
    cType <- base::sapply(strsplit(colnames(Profiles$expr),"\\|"),head,1)
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
  sel <- reshape2::melt(diag(ncol(Profiles$sets)))$value!=1
  ix <- base::sample(nrow(Profiles$sets),min(selN,nrow(Profiles$sets)))
  allG <- unique(unlist(Profiles$selG))
  for(i in ix){
    #if(which(ix==i)%%50==0) cat(which(ix==i),"in",selN,"\n")
    one <- sort(Profiles$sets[i,])
    res <- cor(Profiles$expr[allG,one],Profiles$expr[allG,one],method=corMethod)
    DIST <- as.matrix(dist(t(Profiles$expr[allG,one])))
    dimnames(res) <- base::sapply(dimnames(res),function(x){return(list(sapply(strsplit(x,"\\|"),head,1)))})
    allCor <- rbind(allCor,cbind(reshape2::melt(res)[sel,],eucDist=reshape2::melt(DIST)[sel,"value"]))
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
    dimnames(res) <- base::sapply(dimnames(res),function(x){return(list(sapply(strsplit(x,"\\|"),head,1)))})
    Data <- rbind(allCor[allCor[,1]==i,],cbind(reshape2::melt(res)[reshape2::melt(lower.tri(res))$value,],eucDist=reshape2::melt(DIST)[reshape2::melt(lower.tri(DIST))$value,"value"]))
    print(ggplot2::ggplot(Data,ggplot2::aes(Var2,value,fill=Var2))+ggplot2::geom_boxplot(notch=T)+
    	  	ggplot2::scale_fill_manual(values=cellCol)+
    	  	ggplot2::scale_y_continuous(trans='log10')+
    	  	ggplot2::ggtitle(paste(main,i,sep=":"))+ggplot2::ylab("Correlation Coefficients")+
    	  	ggplot2::theme_classic()+
    	  	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)))
    print(ggplot2::ggplot(Data,ggplot2::aes(Var2,eucDist,fill=Var2))+ggplot2::geom_boxplot(notch=T)+
    	  	ggplot2::scale_fill_manual(values=cellCol)+
    	  	ggplot2::scale_y_continuous(trans='log10')+
    	  	ggplot2::ggtitle(paste(main,i,sep=":"))+ggplot2::ylab("Euclidean Distance")+
    	  	ggplot2::theme_classic()+
    	  	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)))
    DS <- setNames(c(DS,popDis(Data)),c(names(DS),i))
  }
  print(ggplot2::ggplot(cbind(as.data.frame(DS),method="NNLS",cType=names(DS)),
  					  ggplot2::aes(x=method,y=DS))+
  	  	ggplot2::geom_boxplot(outlier.shape = NA)+
  	  	ggplot2::geom_point(ggplot2::aes(color=cType),position=ggplot2::position_jitterdodge(jitter.width=0.2),size=2)+
  	  	ggplot2::scale_color_manual(values=cellCol)+
  	  	ggplot2::ylab("DS")+
  	  	ggplot2::theme_classic())
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
