##############################
## plotPCA.R
##
## marker: a list of two, first one with a column name in pheno indicated the color;
##                        second ne with a column name in pheno indicated the shape
#############################
if(!suppressWarnings(suppressMessages(require(ggfortify)))) install.packages("ggfortify",repos="https://cloud.r-project.org/")
if(!suppressWarnings(suppressMessages(require(ggrepel)))) install.packages("ggrepel",repos="https://cloud.r-project.org/")
if(!require(ggfortify)||!require(ggrepel)) stop("ggfortify or ggrepel cannot be loaded!")


printPCA <- function(X,pheno,signs,main){
  print(plotPCA(X,pheno[colnames(X),],signs)+
          ggtitle(paste(main,nrow(X),"features"))+
          theme(legend.key.size=unit(0.5,"cm")))
}


plotPCA <- function(eData,pheno,marker){
  if(length(marker)==1){
    p <- autoplot(prcomp(t(eData)),data=pheno,colour=names(marker)) + 
      scale_color_manual(values=marker[[1]])+
      theme_classic()
  }else if(length(marker)==2){
    p<- autoplot(prcomp(t(eData)),data=pheno,colour=names(marker)[1],shape=names(marker)[2]) + 
      scale_color_manual(values=marker[[1]])+
      scale_shape_manual(values=marker[[2]])+ 
      theme_classic()
  }else{
    stop("Unknown marker information")
  }
  return(p)
}
projPCA <- function(eData,projD,pheno,marker){
  pca <- prcomp(t(eData))
  prjPCA <- cbind(predict(pca,t(projD)),pheno)
  if(length(marker)==1){
    p <- autoplot(pca$x,colour="gray90") +
      geom_point(data=prjPCA,mapping=aes_string(x="PC1",y="PC2",colour=names(marker)))+
      scale_color_manual(values=marker[[1]])+
      theme_classic()
  }else if(length(marker)==2){
    p<- ggplot() + #,colour="gray90",inherit.aes = F
      geom_point(data=pca$x,mapping=aes(x=PC1,y=PC2),colour="gray90")+
      geom_point(data=prjPCA,mapping=aes_string(x="PC1",y="PC2",colour=names(marker)[1],shape=names(marker)[2]))+
      scale_color_manual(values=marker[[1]])+
      scale_shape_manual(values=marker[[2]])+ 
      theme_classic()
  }else{
    stop("Unknown marker information")
  }
  return(p)
}
projectPCA <- function(eData,pheno,marker,projD){
  pca <- prcomp(t(eData))
  prjPCA <- predict(pca,t(projD))
  pca.var <- summary(pca)$importance
  if(length(marker)==1){
    p<- ggplot() + labs(x=paste("PC1: ",round(pca.var[2,"PC1"]*100,2),"%",sep=""),
                        y=paste("PC2: ",round(pca.var[2,"PC2"]*100,2),"%",sep="")) +
      geom_point(data=cbind(pca$x,pheno),mapping=aes_string(x="PC1",y="PC2",colour=names(marker)))+
      scale_color_manual(values=marker[[1]])+
      geom_point(data=prjPCA,mapping=aes(x=PC1,y=PC2),color="#000000")+
      geom_text_repel(data=prjPCA,aes(x=PC1,y=PC2,label=rownames(prjPCA)),size=2)+
      theme_classic()
  }else if(length(marker)==2){
    p<- ggplot() + #,colour="gray90",inherit.aes = F
      geom_point(data=cbind(pca$x,pheno),mapping=aes_string(x="PC1",y="PC2",colour=names(marker)[1],shape=names(marker)[2]))+
      scale_color_manual(values=marker[[1]])+
      scale_shape_manual(values=marker[[2]])+ 
      geom_point(data=prjPCA,mapping=aes(x=PC1,y=PC2),color="#000000")+
      geom_text_repel(data=prjPCA,aes(x=PC1,y=PC2,label=rownames(prjPCA)),size=2)+
      theme_classic()
  }else{
    stop("Unknown marker information")
  }
  return(p)
}


