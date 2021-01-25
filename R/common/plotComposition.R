##############################
## plotComposition.R
##
## comp and pV [cell type,sample index]
##############################

plotComposition <- function(comp,cellCol,pV=NULL,
                            alpha=c("8F","BF","DF","FF"),
                            pSig=NULL,#c(0.05,0.01,0.001,1e-5)
                            scale=T,cutoff.comp=0.05,cutoff.p=0.05,
                            shading=F,ggplotFun=NULL,
                            pReturn=F,rateReturn=F){
  if(scale) comp <- base::apply(comp,2,function(x){return(x/sum(x))})
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
    comp <- base::sapply(colnames(ix),function(sID,X){
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
  pCol <- base::apply(pV,2,function(x){
    col <- cellCol[names(x)]
    logP <- -log10(1e-300+x)
    steps <- 0:length(alpha)*(max(logP)+1)/length(alpha)
    col <- paste(col,alpha[sapply(logP,function(x){return(which((steps-x)>0)[1]-1)})],sep="")
    return(col)
  })
  pTxt <- base::apply(pV,2,function(x){
    txt <- setNames(rep("",length(x)),names(x))
    for(i in 1:length(pSig)){
      txt[x<pSig[i]] <- paste(rep("*",i),collapse="")
    }
    return(txt)
  })
  pAst <- base::sapply(1:length(pSig),function(x){return(paste(rep("*",x),collapse=""))})
  pAst <- setNames(paste(pAst,"<",format(pSig,digits=1,scientific=T),sep=""),pAst)
  ## plotting -------
  D <- cbind(reshape2::melt(comp),COL=reshape2::melt(pCol)$value,pTxt=reshape2::melt(pTxt)$value)
  p <- ggplot2::ggplot(data=D,ggplot2::aes(x=as.character(Var2),y=value,fill=Var1))+
  	ggplot2::geom_bar(stat="identity")+
  	ggplot2::scale_fill_manual(values=cellCol,name="Cell Type") +
  	ggplot2::ylab("Composition")+
  	ggplot2::theme_classic()+
  	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
          axis.title.x=ggplot2::element_blank()#,axis.title.y=element_blank()
          )
  if(sum(nchar(pTxt)>0)>0) p <- p+
  	ggplot2::geom_text(ggplot2::aes(label=pTxt),position = ggplot2::position_stack(vjust = 0.5)) +
  	ggplot2::geom_point(ggplot2::aes(x=x,y=y,shape=pValue),
               data=data.frame(x=rep(1,length(pAst)),
                               y=rep(0,length(pAst)),
                               pValue=factor(pAst,levels=pAst),
                               Var1=rep(names(cellCol)[1],length(pAst))))+
  	ggplot2::scale_shape_manual(values=setNames(rep(NA,length(pAst)+length(cellCol)),c(pAst,names(cellCol))))+
  	ggplot2::guides(shape=ggplot2::guide_legend(title.hjust = 1,order=0),
  					fill=ggplot2::guide_legend(override.aes=list(shape=NA),order=1))
  if(is.null(ggplotFun)) print(p)
  else print(p+ggplotFun)
  for(selC in rownames(comp)){
    if(selC=="undetermined" || max(comp[selC,])<0.05) next
    fig <- ggplot2::ggplot(data=D[D$Var1==selC,],ggplot2::aes(x=as.character(Var2),y=value,fill=Var1))+
    	ggplot2::geom_bar(stat="identity")+
    	ggplot2::scale_fill_manual(values=cellCol[selC],name="Cell Type") +
    	ggplot2::ylab("Composition")+ggplot2::ggtitle(selC)+ggplot2::ylim(0,1)+
    	ggplot2::theme_classic()+
    	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,face="bold",size=10),
            axis.title.x=ggplot2::element_blank()#,axis.title.y=element_blank()
      )
    if(sum(nchar(pTxt)>0)>0) fig <- fig+
    		ggplot2::geom_text(ggplot2::aes(label=pTxt),position = ggplot2::position_stack(vjust = 0.5)) +
    		ggplot2::geom_point(ggplot2::aes(x=x,y=y,shape=pValue),
                   data=data.frame(x=rep(1,length(pAst)),
                                   y=rep(0,length(pAst)),
                                   pValue=factor(pAst,levels=pAst),
                                   Var1=rep(selC,length(pAst))))+
    		ggplot2::scale_shape_manual(values=setNames(rep(NA,length(pAst)+1),c(pAst,selC)))+
    		ggplot2::guides(shape=ggplot2::guide_legend(title.hjust = 1,order=0),
    						fill=ggplot2::guide_legend(override.aes=list(shape=NA),order=1))
    print(fig)
  }
  
  if(shading){
    p <- ggplot2::ggplot(data=D,ggplot2::aes(x=Var2,y=value,fill=COL))+
    	ggplot2::geom_bar(stat="identity")+
    	ggplot2::scale_fill_identity() +
    	ggplot2::theme_classic()+
    	ggplot2::theme(legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            axis.title.x=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank())
    if(is.null(ggplotFun)) print(p)
    else print(p+ggplotFun)
  }
  if(pReturn) return(p+ggplotFun)
  if(rateReturn) return(comp[rownames(comp)!="undetermined",])
}
plotCorComposition <- function(predictR,mixR,cellCol){
  predictR <- base::apply(predictR,2,function(x){return(x/sum(x))})
  coeff <- diag(cor(predictR,mixR))
  #print(coeff)
  coeff[!is.finite(coeff)] <- 0
  ## plotting -----
  D <- cbind(reshape2::melt(mixR),predict=reshape2::melt(predictR)$value)
  colnames(D)[1:3] <- c("cellType","sample ID","expect")
  print(ggplot2::ggplot(D,ggplot2::aes(x=expect,y=predict,colour=cellType))+
  	  	ggplot2::scale_color_manual(values = cellCol) +
  	  	ggplot2::ggtitle("Overall") +
  	  	ggplot2::geom_point()+ggplot2::theme_classic())
  for(i in colnames(predictR)){
    D <- data.frame(expect=mixR[,i],
                    predict=predictR[,i],
                    cellType=rownames(predictR))
    print(ggplot2::ggplot(D,ggplot2::aes(x=expect,y=predict,colour=cellType))+
    	  	ggplot2::scale_color_manual(values = cellCol) +
    	  	ggplot2::ggtitle(paste(i,round(cor(mixR[,i],predictR[,i]),2),sep=":")) +
    	  	ggplot2::geom_point()+ggplot2::theme_classic())
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
