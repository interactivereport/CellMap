require(ggpubr)
require(ggplot2)
require(reshape2)
require(scales)
require(ggfortify)

rm(list=ls())
graphics.off()
closeAllConnections()

Profiles <- c("CNS6","Neuron3")#
PCA <- readRDS("S1.2.rds")
techq <- c('DropSeq','SMARTseq','Bulk','fluidigm C1','10X')
techCol <- c(setNames(dscale(techq, hue_pal(l = 75)),techq),
			 setNames(dscale(factor(1:2), hue_pal(l = 75)),c('snRNA','scRNA')))

plotPCA <- function(onePCA,pheno,marker){
	if(length(marker)==1){
		p <- autoplot(onePCA,data=pheno,colour=names(marker)) + 
			scale_color_manual(values=marker[[1]])
	}else if(length(marker)==2){
		p<- autoplot(onePCA,data=pheno,colour=names(marker)[1],shape=names(marker)[2]) + 
			scale_color_manual(values=marker[[1]])+
			scale_shape_manual(values=marker[[2]])
	}
	p <- p+theme_light()+
		theme(axis.text = element_text(size=10),
			  axis.title=element_text(size=11),
			  aspect.ratio=1,
			  panel.grid.minor=element_blank(),
			  panel.grid.major=element_line(colour="#d9d9d9",linetype="dashed"),
			  plot.margin=unit(c(0.1,0.1,0.1,0.1),"inches"),
			  legend.position = "none")
	return(p)
}
plotDist <- function(Data,col,facetLabel){
	p <- ggplot(Data,aes(x=expr,color=sID))+
		facet_wrap(~cellType,ncol=3,labeller=facetLabel)+
		geom_density(geom='raster')+
		scale_color_manual(values=col)+
		scale_x_continuous(trans = "identity",limits=c(1,15))+
		scale_y_continuous(limits=c(0,0.4))+
		ylab("Density") +
		theme_light()+
		theme(legend.position = "none",
			  axis.text = element_text(size=10),
			  axis.title=element_text(size=10),
			  aspect.ratio=1,
			  strip.text = element_text(size = 7,color="black"),
			  strip.background = element_rect(fill="white",colour="gray"),
			  panel.grid.minor=element_blank(),
			  panel.grid.major=element_line(colour="#d9d9d9",linetype="dashed"),
			  plot.margin=unit(c(0.1,0.1,0.1,0.1),"inches"))
	return(p)
}
scaleLabel <- function(x){
	if(max(x,na.rm=T)<50) return(x)
	a <- format(x,digits=2,scientific=T)
	a[is.na(x)] <- NA
	return(a)
}


wR <- c(1,1,1)
hR <- c(1,0.5)
for(one in Profiles){
	message("=============== ",one," ===============")
	D <- readRDS(paste0("../../profiles/",one,".rds"))
	pseudoID <- rownames(PCA[[one]]$rawPCA$x)
	pheno <- data.frame(row.names=pseudoID,
						cell_type=sapply(strsplit(pseudoID,"\\|"),head,1),
						dataset=sapply(strsplit(pseudoID,"\\|"),"[[",2),
						technology=D$para$tech[sapply(strsplit(pseudoID,"\\|"),"[[",2),1])
	figs <- list()
	## PCA plots -----
	message("plotting PCA ...")
	figs[[length(figs)+1]] <- plotPCA(PCA[[one]]$rawPCA,pheno,
									  list(technology=techCol[levels(pheno$technology)]))
	figs[[length(figs)+1]] <- plotPCA(PCA[[one]]$rawPCA,pheno,
									  list(cell_type=D$para$cellCol[levels(pheno$cell_type)],
									  	   dataset=setNames(1:nlevels(pheno$dataset),levels(pheno$dataset))))
	
	figs[[length(figs)+1]] <- plotPCA(PCA[[one]]$fullPCA,pheno,
									  list(cell_type=D$para$cellCol[levels(pheno$cell_type)],
									  	 dataset=setNames(1:nlevels(pheno$dataset),levels(pheno$dataset))))

	figLabel <- LETTERS[1:length(figs)]
	## obtain the legend -----
	p <- figs[[length(figs)]]+
		geom_col(aes(x,y,fill=technology),
				 data.frame(x=1:nlevels(pheno$technology),
				 		   y=1:nlevels(pheno$technology),
				 		   technology=levels(pheno$technology)))+
		scale_fill_manual(values=techCol[levels(pheno$technology)])+
		theme(legend.position='right',
			  legend.text=element_text(size=8,margin=margin(t=0.5)),
			  legend.title= element_text(size=9),
			  legend.key.size=unit(3,'mm'),
			  legend.spacing = unit(-2,'mm'))+
		guides(fill=guide_legend(ncol=2,order=1),
			   shape=guide_legend(ncol=2,order=2),
			   color=guide_legend(ncol=2,order=3))
	
	figs[[length(figs)+1]] <- as_ggplot(get_legend(p))+
		theme(plot.margin=unit(c(0.2,0.2,0.1,0.2),"inches"))
	figLabel <- c(figLabel,"")

	## Distribution ----
	message("plotting distribution ...")
	X <- melt(as.matrix(D$expr[unique(unlist(D$selG)),]),
			  varnames=c('genes','sID'),value.name='expr')
	X[['cellType']] <- sapply(strsplit(as.character(X$sID),"\\|"),head,1)
	sCol <- setNames(D$para$cellCol[sapply(strsplit(levels(X$sID),"\\|"),head,1)],
					 levels(X$sID))
	
	figs[[length(figs)+1]] <- plotDist(X,sCol,'label_value')+
		xlab(paste0("log2(UMI+1) (",nlevels(X$genes),"genes)"))
	figLabel <- c(figLabel,LETTERS[length(figs)-1])
	for(i in names(D$selG)){
		X <- X[X$cellType!=i | X$genes%in%D$selG[[i]],]
	}
	cellTypeN <- sapply(D$selG,length)
	figs[[length(figs)+1]] <- plotDist(X,sCol,
									   labeller(cellType=setNames(paste(substr(names(cellTypeN),1,3),
									   								 cellTypeN,sep=".:"),
									   						   names(cellTypeN))))+
		xlab("log2(UMI+1)")
	figLabel <- c(figLabel,LETTERS[length(figs)-1])

	## saving -----
	message("saving ...")
	pdf(paste0("S",which(Profiles==one),".pdf"),width=8.5,height=11)
	print(ggarrange(plotlist=figs,labels=figLabel,
					#widths=wR,
					heights=c(1,hR[which(Profiles==one)],1.8+(1-hR[which(Profiles==one)])),
					ncol=3,nrow=3)+
		  	theme(plot.margin=unit(c(0.4,0.1,0.1,0.1),"inches")))

	a <- dev.off()
}



