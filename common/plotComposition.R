##############################
## plotComposition.R
##
## comp and pV [cell type,sample index]
##############################
loadPlotComposition <- function(){
  if(!require(reshape2,quietly=T, warn.conflicts=F)){
    install.packages("reshape2",repos="https://cloud.r-project.org/")
    if(!require(reshape2)) stop("Cannot install package reshape2!")
  }
  if(!require(ggplot2,quietly=T, warn.conflicts=F)){
    install.packages("ggplot2",repos="https://cloud.r-project.org/")
    if(!require(ggplot2)) stop("Cannot install package ggplot2!")
  }
}
suppressWarnings(suppressMessages(loadPlotComposition()))

plotComposition <- function(comp,cellCol,pV=NULL,
                            alpha=c("8F","BF","DF","FF"),
                            pSig=NULL,#c(0.05,0.01,0.001,1e-5)
                            scale=T,cutoff.comp=0.05,cutoff.p=0.05,
                            shading=F,ggplotFun=NULL,
                            pReturn=F,rateReturn=F){
  if(scale) comp <- apply(comp,2,function(x){return(x/sum(x))})
  comp <- comp[order(rownames(comp)),]
  if(is.null(pV)){
    pV <- comp
    pV[,] <- 0.9
    if(is.null(pSig))
      pSig <- c(0.05,0.01,0.001,1e-5)
  }else{
    pV <- pV[rownames(comp),]
    ## addding undetermined ------
    cellCol <- c(cellCol,undetermined="#BBBBBB")
    pV[comp<cutoff.comp] <- 1
    ix <- pV > cutoff.p
    comp <- sapply(colnames(ix),function(sID,X){
      x <- setNames(c(X[,sID],sum(X[ix[,sID],sID])),c(rownames(X),"undetermined"))
      x[c(ix[,sID],F)] <- 0
      return(x)
    },comp)
    pV <- rbind(pV,undetermined=rep(1,ncol(pV)))
    ## dynamic determin the p-value asterisk
    if(is.null(pSig))
      pSig <- unique(as.numeric(format(c(cutoff.p,summary(pV[pV<cutoff.p])[c(5,3,2)]),digits=1,scientific=T)))
    if(length(pSig)<4) pSig <- c(0.05,0.01,0.001,1e-5)
    if(pSig[4]==0) pSig[4] <- as.numeric(format(min(pV[pV>0]),digits=1,scientific=T))
  }
  ## shade the p-value or * them -----
  pCol <- apply(pV,2,function(x){
    col <- cellCol[names(x)]
    logP <- -log10(1e-300+x)
    steps <- 0:length(alpha)*(max(logP)+1)/length(alpha)
    col <- paste(col,alpha[sapply(logP,function(x){return(which((steps-x)>0)[1]-1)})],sep="")
    return(col)
  })
  pTxt <- apply(pV,2,function(x){
    txt <- setNames(rep("",length(x)),names(x))
    for(i in 1:length(pSig)){
      txt[x<pSig[i]] <- paste(rep("*",i),collapse="")
    }
    return(txt)
  })
  pAst <- sapply(1:length(pSig),function(x){return(paste(rep("*",x),collapse=""))})
  pAst <- setNames(paste(pAst,"<",format(pSig,digits=1,scientific=T),sep=""),pAst)
  ## plotting -------
  D <- cbind(melt(comp),COL=melt(pCol)$value,pTxt=melt(pTxt)$value)
  p <- ggplot(data=D,aes(x=as.character(Var2),y=value,fill=Var1))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=cellCol,name="Cell Type") +
    ylab("Composition")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x=element_blank()#,axis.title.y=element_blank()
          )
  if(sum(nchar(pTxt)>0)>0) p <- p+
    geom_text(aes(label=pTxt),position = position_stack(vjust = 0.5)) +
    geom_point(aes(x=x,y=y,shape=pValue),
               data=data.frame(x=rep(1,length(pAst)),
                               y=rep(0,length(pAst)),
                               pValue=factor(pAst,levels=pAst),
                               Var1=rep(names(cellCol)[1],length(pAst))))+
    scale_shape_manual(values=setNames(rep(NA,length(pAst)+length(cellCol)),c(pAst,names(cellCol))))+
    guides(shape=guide_legend(title.hjust = 1,order=0),fill=guide_legend(override.aes=list(shape=NA),order=1))
  if(is.null(ggplotFun)) print(p)
  else print(p+ggplotFun)
  for(selC in rownames(comp)){
    if(selC=="undetermined" || max(comp[selC,])<0.05) next
    fig <- ggplot(data=D[D$Var1==selC,],aes(x=as.character(Var2),y=value,fill=Var1))+
      geom_bar(stat="identity")+
      scale_fill_manual(values=cellCol[selC],name="Cell Type") +
      ylab("Composition")+ggtitle(selC)+ylim(0,1)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,face="bold",size=10),
            axis.title.x=element_blank()#,axis.title.y=element_blank()
      )
    if(sum(nchar(pTxt)>0)>0) fig <- fig+
        geom_text(aes(label=pTxt),position = position_stack(vjust = 0.5)) +
        geom_point(aes(x=x,y=y,shape=pValue),
                   data=data.frame(x=rep(1,length(pAst)),
                                   y=rep(0,length(pAst)),
                                   pValue=factor(pAst,levels=pAst),
                                   Var1=rep(selC,length(pAst))))+
        scale_shape_manual(values=setNames(rep(NA,length(pAst)+1),c(pAst,selC)))+
        guides(shape=guide_legend(title.hjust = 1,order=0),fill=guide_legend(override.aes=list(shape=NA),order=1))
    print(fig)
  }
  
  if(shading){
    p <- ggplot(data=D,aes(x=Var2,y=value,fill=COL))+
      geom_bar(stat="identity")+
      scale_fill_identity() +
      theme_classic()+
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    if(is.null(ggplotFun)) print(p)
    else print(p+ggplotFun)
  }
  if(pReturn) return(p+ggplotFun)
  if(rateReturn) return(comp[rownames(comp)!="undetermined",])
}
plotCorComposition <- function(predictR,mixR,cellCol){
  predictR <- apply(predictR,2,function(x){return(x/sum(x))})
  coeff <- diag(cor(predictR,mixR))
  #print(coeff)
  coeff[!is.finite(coeff)] <- 0
  ## plotting -----
  D <- cbind(melt(mixR),predict=melt(predictR)$value)
  colnames(D)[1:3] <- c("cellType","sample ID","expect")
  print(ggplot(D,aes(x=expect,y=predict,colour=cellType))+
          scale_color_manual(values = cellCol) +
          ggtitle("Overall") +
          geom_point()+theme_classic())
  for(i in colnames(predictR)){
    D <- data.frame(expect=mixR[,i],
                    predict=predictR[,i],
                    cellType=rownames(predictR))
    print(ggplot(D,aes(x=expect,y=predict,colour=cellType))+
            scale_color_manual(values = cellCol) +
            ggtitle(paste(i,round(cor(mixR[,i],predictR[,i]),2),sep=":")) +
            geom_point()+theme_classic())
  }
  return(coeff)
}
getRMSE <- function(predictR,mixR){
  if(is.null(ncol(predictR))){
    return(sqrt(mean((predictR-mixR)^2)))
  }
  rmse <- apply(mixR-predictR[rownames(mixR),colnames(mixR)],2,function(x){return(sqrt(mean(x^2)))})
  print(rmse)
  return(rmse)
}
getCOR <- function(predictR,mixR,method="pearson"){
  if(is.null(ncol(predictR))){
    return(cor(predictR,mixR,method=method))
  }
  return(diag(cor(mixR,predictR[rownames(mixR),colnames(mixR)],method=method)))
  
}
## weighted RMSE, penalize the false positive population
getWRMSE <- function(predictR,mixR,penalty=2){
  w <- mixR
  w[mixR>0] <- 1
  w[mixR==0] <- penalty
  return(apply(w*(mixR-predictR[rownames(mixR),colnames(mixR)]),2,function(x){return(sqrt(mean(x^2)))}))
}
