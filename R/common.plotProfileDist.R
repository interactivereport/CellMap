#########################
## plotProfileDist.R
##
########################

plotProfileDist <- function(Exp,cellCOL,selG=NULL,...){
  strTitle <- paste("Profile distribution on",nrow(Exp),"genes")
  Data <- reshape2::melt(as.matrix(Exp))
  plotScale <- "identity"
  plotLimits <- c(1,20)
  if(max(Exp)>50){
    plotScale <- "log2"
    plotLimits <- c(2,round(max(Exp)/1000)*1000)
  }
  Data <- cbind(Data,cellType=sapply(strsplit(as.character(Data$Var2),"\\|"),head,1))
  sCol <- setNames(cellCOL[sapply(strsplit(colnames(Exp),"\\|"),head,1)],colnames(Exp))
  figures <- list()
  figures[[length(figures)+1]] <- ggplot2::ggplot(Data,ggplot2::aes(x=value,color=Var2))+
  	ggplot2::facet_wrap(~cellType,ncol=3)+
  	ggplot2::geom_density()+
  	ggplot2::scale_color_manual(values=sCol)+
  	ggplot2::scale_x_continuous(trans = plotScale,limits=plotLimits)+
  	ggplot2::scale_y_continuous(limits=c(0,0.5))+
  	ggplot2::xlab("Profile Expression") +
  	ggplot2::ggtitle(strTitle) + 
  	ggplot2::theme_classic()+
  	ggplot2::theme(legend.position = "none",aspect.ratio=1,
          panel.grid.major.x=ggplot2::element_line(colour="#d9d9d9",linetype="dashed"))
  figures[[length(figures)+1]] <- ggplot2::ggplot(Data,ggplot2::aes(x=value,color=cellType))+
  	ggplot2::geom_density()+
  	ggplot2::scale_color_manual(values=cellCOL)+
  	ggplot2::scale_x_continuous(trans = plotScale,limits=plotLimits)+
  	ggplot2::scale_y_continuous(limits=c(0,0.5))+
  	ggplot2::xlab("Profile Expression") +
  	ggplot2::ggtitle(strTitle) + 
  	ggplot2::theme_classic()+
  	ggplot2::theme(legend.position = "none",aspect.ratio=1,
          panel.grid.major.x=ggplot2::element_line(colour="#d9d9d9",linetype="dashed"))
  ## only cell type specific gene
  if(!is.null(selG)){
    annoGene <- matrix("Others",nrow=nrow(Exp),ncol=length(selG),dimnames=list(rownames(Exp),names(selG)))
    annoGeneCol <- list()
    for(i in names(selG)){
      figures[[length(figures)+1]] <- ggplot2::ggplot(Data[Data$Var1%in%selG[[i]],],ggplot2::aes(x=value,color=Var2))+
      	ggplot2::facet_wrap(~cellType,ncol=3)+
      	ggplot2::geom_density()+
      	ggplot2::scale_color_manual(values=sCol)+
      	ggplot2::scale_x_continuous(trans = plotScale,labels=scaleLabel,limits=plotLimits)+#
      	ggplot2::scale_y_continuous(limits=c(0,0.5))+
      	ggplot2::xlab("Profile Expression") +
      	ggplot2::ggtitle(paste(i,": ",length(selG[[i]])," genes",sep="")) + 
      	ggplot2::theme_classic()+
      	ggplot2::theme(legend.position = "none",
              axis.text.x = ggplot2::element_text(angle = 90),
              aspect.ratio=1,
              panel.grid.major.x=ggplot2::element_line(colour="#d9d9d9",linetype="dashed"))
      annoGene[selG[[i]],i] <- i
      annoGeneCol[[i]] <- setNames(c("#FFFFFF",cellCOL[i]),c("Others",i))
        
    }
    cat("Before cell type:",nrow(Data),"\n")
    for(i in names(selG)){
      Data <- Data[!(Data$cellType==i&!Data$Var1%in%selG[[i]]),]
    }
    cat("After cell type:",nrow(Data),"\n")
    cellTypeGN <- sapply(selG,function(x){return(sum(x%in%Data$Var1))})
    figures[[length(figures)+1]] <- ggplot2::ggplot(Data,ggplot2::aes(x=value,color=Var2))+
    	ggplot2::facet_wrap(~cellType,ncol=3)+
    	ggplot2::geom_density()+
    	ggplot2::scale_color_manual(values=sCol)+
    	ggplot2::scale_x_continuous(trans = plotScale,labels=scaleLabel,limits=plotLimits)+#
    	ggplot2::scale_y_continuous(limits=c(0,0.5))+
    	ggplot2::xlab("Profile Expression") +
    	ggplot2::ggtitle("Signature genes for each cell type") + 
    	ggplot2::theme_classic()+
    	ggplot2::theme(legend.position = "none",
            aspect.ratio=1,
            axis.text.x = ggplot2::element_text(angle = 90),
            panel.grid.major.x=ggplot2::element_line(colour="#d9d9d9",linetype="dashed"))
    figures[[length(figures)+1]] <- ggplot2::ggplot(Data,ggplot2::aes(x=value,color=cellType))+
    	ggplot2::geom_density()+
    	ggplot2::scale_color_manual(values=cellCOL)+
    	ggplot2::scale_x_continuous(trans = plotScale,labels=scaleLabel,limits=plotLimits)+#
    	ggplot2::scale_y_continuous(limits=c(0,0.5))+
    	ggplot2::xlab("Profile Expression") +
    	ggplot2::ggtitle("Signature genes for each cell type") + 
    	ggplot2::theme_classic()+
    	ggplot2::theme(legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90),
            aspect.ratio=1)
    ## heatmap
    heatCol <- c("#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
                 "#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5","#F8EFEA","#FBE0D1",
                 "#FAC9B1","#F5AD8C","#E88B6E","#D96752","#C6413E","#B31B2C",
                 "#8E0C25","#67001F")
    steps <- c(0,stats::quantile(Exp[Exp>0],probs=round((1:length(heatCol))/length(heatCol),2)))
    Exp <- Exp[,order(colnames(Exp))]
    pheatmap::pheatmap(t(Exp),color=heatCol,breaks=steps,
             cluster_rows=F,cluster_cols=F,#clustering_method="ward.D",
             show_rownames=F,show_colnames=F,#main=i,
             annotation_row=data.frame(row.names=colnames(Exp),cellType=base::sapply(strsplit(colnames(Exp),"\\|"),head,1)),
             annotation_col = as.data.frame(annoGene),
             annotation_colors=c(list(cellType=cellCOL),annoGeneCol),
             ...)
    
  }
  rm(Data)
  return(figures)
}
plotBulkFeature <- function(Bulk,Exp,cellCOL,selG){
  Exp <- Exp[,order(colnames(Exp))]
  D <- cbind(Exp,Bulk)
  cType <- c(base::sapply(strsplit(colnames(Exp),"\\|"),head,1),
             rep("Other",ncol(Bulk)))
  gType <- rep("Other",nrow(D))
  for(i in names(selG)) gType[rownames(D)%in%selG[[i]]] <- i
  
  sID <- c(rep("",ncol(Exp)),colnames(Bulk))
  heatCol <- c("#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
               "#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5","#F8EFEA","#FBE0D1",
               "#FAC9B1","#F5AD8C","#E88B6E","#D96752","#C6413E","#B31B2C",
               "#8E0C25","#67001F")
  steps <- c(0,stats::quantile(D[D>0],probs=round((1:length(heatCol))/length(heatCol),2)))
  pheatmap::pheatmap(t(D),color=heatCol,breaks=steps,
           cluster_rows=F,cluster_cols=F,#clustering_method="ward.D",
           show_colnames=F,#main=i,show_rownames=F,
           fontsize_row=0.5,labels_row=sID,
           annotation_row=data.frame(row.names=colnames(D),cellType=cType),
           annotation_col = data.frame(row.names=rownames(D),geneType=gType),
           annotation_colors=list(cellType=c(cellCOL,Other="#FFFFFF"),geneType=cellCOL))
  pheatmap::pheatmap(t(Bulk),color=heatCol,breaks=steps,
           cluster_rows=F,cluster_cols=F,#clustering_method="ward.D",
           show_colnames=F,#main=i,show_rownames=F,
           fontsize_row=8,labels_row=colnames(Bulk),
           annotation_col = data.frame(row.names=rownames(D),geneType=gType),
           annotation_colors=list(geneType=cellCOL))#cellType=c(cellCOL,Other="#FFFFFF"),
}

scaleLabel <- function(x){
  if(max(x,na.rm=T)<50) return(x)
  a <- format(x,digits=2,scientific=T)
  a[is.na(x)] <- NA
  return(a)
}
