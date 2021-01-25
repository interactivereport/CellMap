######################
## plotEval.R
##
######################

plotEval <- function(estR,expR,cellCol,sID=F,pCH=NULL,rmseMax=0.5){
  if(sum(dim(estR)!=dim(expR))>0) stop("estR and expR are not matching dimension in plotEval")
  if(sum(rownames(estR)!=rownames(expR))>0 || sum(colnames(estR)!=colnames(expR))>0)
    stop("estR and expR are not matching names in plotEval")
  D <- reshape2::melt(as.matrix(estR))
  D <- setNames(reshape2::melt(as.matrix(estR)),c('cellType','sID','Estimated'))
  if(!is.null(pCH)) D <- cbind(D,PCH=pCH[D$sID])

  D$sID <- as.factor(D$sID)
  D <- cbind(D,Expected=reshape2::melt(as.matrix(expR))[,'value'])
  plotEvalOne(D,'cellType',cellCol)
  rmse <- plotRMSE(D,cellCol,yMax=rmseMax)
  #plotEvalSample(D,cellCol,1:9,ggtitle("top 9"))
  #plotEvalSample(D,cellCol,10:18,ggtitle("top 10~8"))
  #plotEvalSample(D,cellCol,19:27,ggtitle("top 19~27"))
  if(sID){
    D$dID <- base::sapply(strsplit(as.character(D$sID),"\\|"),head,1)
    plotEvalOne(D,'dID')
    for(i in unique(D$dID))
      plotEvalOne(D[D$dID==i,],'cellType',cellCol,ggplot2::ggtitle(i))
  }
  return(rmse)
}
plotEvalOne <- function(D,colName,col=NULL,ggFun=NULL,facet=T,Plot=T){
  Dstat <- merge(plyr::ddply(D,colName,dplyr::summarise,
                 Expected=max(Expected),
                 Estimated=max(Estimated)),
  			   plyr::ddply(D,colName,function(x)return(data.frame(rho=round(cor(x$Expected,x$Estimated),2)))),
                 by=colName,all=T)
  p <- ggplot2::ggplot(D,ggplot2::aes_string('Expected','Estimated',color=colName))
  if(colName=='cellType' && "PCH"%in%colnames(D)) p <- p+ ggplot2::geom_point(ggplot2::aes(shape=PCH),size=2)
  else p <- p + ggplot2::geom_point(size=2)
  p <- p+
  	ggplot2::geom_smooth(method=lm,fill=NA)+
  	ggplot2::geom_abline(intercept=0, slope=1, linetype='dashed', color='gray')+
  	ggplot2::theme_classic(base_size=15)+ggplot2::theme(aspect.ratio=1)
  if(!is.null(col)) p <- p + ggplot2::scale_color_manual(values=col)
  if(!is.null(ggFun)) p <- p + ggFun
  if(Plot){
    print(p+ggrepel::geom_label_repel(ggplot2::aes(label=rho),data=Dstat,
                             box.padding=0.1,size=2,
                             segment.color = 'grey50'))
    if(facet){
      p <- p+ggplot2::xlim(0,1)+ggplot2::ylim(0,1)+
      	ggplot2::facet_wrap(as.formula(paste("~",colName)))+
      	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,hjust=1))
      print(p)
      Dstat[,"Expected"] <- 0 
      Dstat[,"Estimated"] <- 1 
      print(p+ggrepel::geom_label_repel(ggplot2::aes(label=rho),data=Dstat,
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

plotRMSE <- function(D,col,ggFun=NULL,yMax=0.5){
  summaD <- plyr::ddply(D,'cellType',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  tmp <- plyr::ddply(D,'sID',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  tmp[,'sID'] <- 'samples'
  colnames(tmp) <- colnames(summaD) <- c("Obv",'RMSE')
  summaD <- rbind(summaD,tmp)
  
  col <- c(col,samples='gray')
  X <- plyr::ddply(summaD,'Obv',dplyr::summarise,
             mu=mean(RMSE),
             sigma=sd(RMSE))
  X[is.na(X)] <- 0
  
  p <- ggplot2::ggplot(X,ggplot2::aes(Obv,mu,fill=Obv))+
  	ggplot2::geom_bar(stat='identity',color='black')+
  	ggplot2::geom_errorbar(ggplot2::aes(ymin=mu,ymax=mu+sigma),width=0.3)+
  	ggplot2::geom_point(ggplot2::aes(Obv,RMSE),data=summaD,position=ggplot2::position_jitterdodge(),pch=21,size=1)+
    #geom_dotplot(aes(Obv,RMSE),data=summaD,binaxis='y',stackdir="center",dotsize=0.1,color="gray50",fill="gray50")+
    #ylim(0,yMax)+
  	ggplot2::coord_cartesian(ylim=c(0, yMax))+
  	ggplot2::scale_fill_manual(values=col)+
  	ggplot2::ylab('RMSE')+ggplot2::xlab("")+
  	ggplot2::theme_classic()+
  	ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=1))
  if(!is.null(ggFun)) p <- p + ggFun
  print(p)
  return(tmp)
}

plotEvalSample <- function(D,col,ix=1:9,ggFun=NULL){
  summaD <- plyr::ddply(D,'sID',function(x)return(data.frame(rmse=sqrt(mean((x$Estimated-x$Expected)^2)))))
  if(nrow(summaD)<min(ix)) return()
  if(nrow(summaD)<max(ix)) ix<- ix[1]:nrow(summaD)
  D1 <- D[as.character(D$sID)%in%summaD[order(summaD$rmse,decreasing=T),'sID'][ix],]
  pch <- setNames(10:(9+length(ix)),unique(D1$sID))
  
  p <- ggplot2::ggplot(D1,ggplot2::aes(Expected,Estimated,color=cellType,fill=cellType,shape=sID))+
  	ggplot2::geom_point()+
  	ggplot2::scale_color_manual(values=col)+
  	ggplot2::scale_shape_manual(values=pch)+
  	ggplot2::geom_abline(intercept=0, slope=1, linetype='dashed', color='gray')+
  	ggplot2::theme_classic()
  if(!is.null(ggFun)) p <- p + ggFun
  print(p)
  print(p+ggplot2::facet_wrap(~sID)+ggplot2::theme(legend.position = "none"))
}




