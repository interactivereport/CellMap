######################
## plotEval.R
##
######################

require(reshape2)
require(ggplot2)
require(ggrepel)
require(plyr)

plotEval <- function(estR,expR,cellCol,sID=F,pCH=NULL){
  if(sum(dim(estR)!=dim(expR))>0) stop("estR and expR are not matching dimension in plotEval")
  if(sum(rownames(estR)!=rownames(expR))>0 || sum(colnames(estR)!=colnames(expR))>0)
    stop("estR and expR are not matching names in plotEval")
  D <- melt(as.matrix(estR))
  D <- setNames(melt(as.matrix(estR)),c('cellType','sID','Estimated'))
  if(!is.null(pCH)) D <- cbind(D,PCH=pCH[D$sID])

  D$sID <- as.factor(D$sID)
  D <- cbind(D,Expected=melt(as.matrix(expR))[,'value'])
  plotEvalOne(D,'cellType',cellCol)
  plotRMSE(D,cellCol)
  plotEvalSample(D,cellCol,1:9,ggtitle("top 9"))
  plotEvalSample(D,cellCol,10:18,ggtitle("top 10~8"))
  plotEvalSample(D,cellCol,19:27,ggtitle("top 19~27"))
  if(sID){
    D$dID <- sapply(strsplit(as.character(D$sID),"\\|"),head,1)
    plotEvalOne(D,'dID')
    for(i in unique(D$dID))
      plotEvalOne(D[D$dID==i,],'cellType',cellCol,ggtitle(i))
  }
}
plotEvalOne <- function(D,colName,col=NULL,ggFun=NULL,facet=T,Plot=T){
  Dstat <- merge(ddply(D,colName,summarise,
                 Expected=max(Expected),
                 Estimated=max(Estimated)),
                 ddply(D,colName,function(x)return(data.frame(rho=round(cor(x$Expected,x$Estimated),2)))),
                 by=colName,all=T)
  p <- ggplot(D,aes_string('Expected','Estimated',color=colName))
  if(colName=='cellType' && "PCH"%in%colnames(D)) p <- p+ geom_point(aes(shape=PCH),size=2)
  else p <- p + geom_point(size=2)
  p <- p+
    geom_smooth(method=lm,fill=NA)+
    geom_abline(intercept=0, slope=1, linetype='dashed', color='gray')+
    theme_classic()+theme(aspect.ratio=1)
  if(!is.null(col)) p <- p + scale_color_manual(values=col)
  if(!is.null(ggFun)) p <- p + ggFun
  if(Plot){
    print(p+geom_label_repel(aes(label=rho),data=Dstat,
                             box.padding=0.1,size=2,
                             segment.color = 'grey50'))
    if(facet){
      p <- p+xlim(0,1)+ylim(0,1)+
        facet_wrap(as.formula(paste("~",colName)))+
        theme(axis.text.x = element_text(angle = 90,hjust=1))
      print(p)
      Dstat[,"Expected"] <- 0 
      Dstat[,"Estimated"] <- 1 
      print(p+geom_label_repel(aes(label=rho),data=Dstat,
                               box.padding=0,label.padding=0.1,size=3,
                               segment.color = 'grey50'))
    }
  }else{
    return(p)
  }
}

test <- function(rmse){
  message(rmse)
  return(0)
}

plotRMSE <- function(D,col,ggFun=NULL){
  summaD <- ddply(D,'cellType',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  tmp <- ddply(D,'sID',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  tmp[,'sID'] <- 'samples'
  colnames(tmp) <- colnames(summaD) <- c("Obv",'RMSE')
  summaD <- rbind(summaD,tmp)
  
  col <- c(col,samples='gray')
  X <- ddply(summaD,'Obv',summarise,
             mu=mean(RMSE),
             sigma=sd(RMSE))
  X[is.na(X)] <- 0
  
  p <- ggplot(X,aes(Obv,mu,fill=Obv))+
    geom_bar(stat='identity',color='black')+
    geom_errorbar(aes(ymin=mu,ymax=mu+sigma),width=0.3)+
    geom_point(aes(Obv,RMSE),data=summaD,position=position_jitterdodge(),pch=21,size=1)+
    #geom_dotplot(aes(Obv,RMSE),data=summaD,binaxis='y',stackdir="center",dotsize=0.1,color="gray50",fill="gray50")+
    ylim(0,0.5)+
    scale_fill_manual(values=col)+
    ylab('RMSE')+xlab("")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,hjust=1))
  if(!is.null(ggFun)) p <- p + ggFun
  print(p)
}

plotEvalSample <- function(D,col,ix=1:9,ggFun=NULL){
  summaD <- ddply(D,'sID',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  if(nrow(summaD)<min(ix)) return()
  if(nrow(summaD)<max(ix)) ix<- ix[1]:nrow(summaD)
  D1 <- D[as.character(D$sID)%in%summaD[order(summaD$rmse,decreasing=T),'sID'][ix],]
  pch <- setNames(10:(9+length(ix)),unique(D1$sID))
  
  p <- ggplot(D1,aes(Expected,Estimated,color=cellType,fill=cellType,shape=sID))+
    geom_point()+
    scale_color_manual(values=col)+
    scale_shape_manual(values=pch)+
    geom_abline(intercept=0, slope=1, linetype='dashed', color='gray')+
    theme_classic()
  if(!is.null(ggFun)) p <- p + ggFun
  print(p)
  print(p+facet_wrap(~sID)+theme(legend.position = "none"))
}




