require(reshape2)
#require(plyr)
require(ggpubr)
require(ggplot2)
require(ggrepel)

rm(list=ls())
graphics.off()
closeAllConnections()

D <- readRDS('Fig4.rds')
ymax <- setNames(c(0.3,0.4,0.5),names(D))
avgData <- F
figs <- list()
bw <- 2.2
wR <- c(0.8,1.2)
for(one in names(D)){
	message(one)
	cellCol <- readRDS(paste0("../../profiles/",one,".rds"))$para$cellCol

	totalX <- NULL
	for(oneT in names(D[[one]])){
		if(sum(dim(D[[one]][[oneT]]$estR)!=dim(D[[one]][[oneT]]$expR))>0 || 
		   sum(rownames(D[[one]][[oneT]]$estR)!=rownames(D[[one]][[oneT]]$expR))>0 ||
		   sum(colnames(D[[one]][[oneT]]$estR)!=colnames(D[[one]][[oneT]]$expR))>0)
			stop(paste("Error non-matching:",oneT,one))
		X <- cbind(melt(D[[one]][[oneT]]$estR,varnames=c('celltype','sID'),value.name='estR'),
				   expR=melt(D[[one]][[oneT]]$expR)$value)
		if(oneT=='pseudo'){
			X <- cbind(X,dID=sapply(strsplit(as.character(X$sID),"\\|"),head,1))
			if(avgData){
				plotX <- plyr::ddply(X,"celltype",function(x)return(data.frame(rmse=sqrt(mean((x$expR-x$estR)^2)))))
				figs[[length(figs)+1]] <- ggplot(plotX,aes(celltype,rmse,fill=celltype))+
					geom_bar(stat="identity")+
					coord_cartesian(ylim=c(0, ymax[one]))+
					ylab('RMSE')+xlab("")+ggtitle(paste(one,"by cell type through all datasets"))+
					scale_fill_manual(values=cellCol)+
					theme_light()+
					theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=12),
						  axis.text.y = element_text(size=12),
						  axis.title=element_text(size=12),
						  legend.position = "none",
						  panel.grid.minor=element_blank(),
						  panel.grid.major.x=element_blank(),
						  plot.margin=unit(c(0.2,0,0.1,0.1),"inches"))
			}
			else{
				plotX <- plyr::ddply(X,c("celltype","dID"),function(x)return(data.frame(rmse=sqrt(mean((x$expR-x$estR)^2)))))
				labelX <- plotX[NULL,]
				for(i in levels(plotX$celltype))labelX <- rbind(labelX,head(plotX[intersect(order(plotX$rmse,decreasing=T),grep(i,plotX$celltype)),],1))
				labelX <- labelX[labelX$rmse>0.05,]
				#plotX$celltype <- factor(plotX$celltype,levels=names(cellCol))
				figs[[length(figs)+1]] <- ggplot(plotX,aes(celltype,rmse,fill=celltype))+
					geom_boxplot(outlier.shape=NA,width=bw*nlevels(plotX$celltype)/20*(sum(wR)-wR[1])/sum(wR))+
					geom_point(position=position_jitterdodge(),pch=21,size=1.5,color="gray")+
					#geom_jitter(alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.2))+
					coord_cartesian(ylim=c(0, ymax[one]))+
					ylab('RMSE')+xlab("")+ggtitle(paste(one,"by cell type through\neach dataset"))+
					scale_fill_manual(values=cellCol)+
					#geom_text_repel(aes(label=dID),data=labelX,size=2)+
					theme_light()+
					theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=12),
						  axis.text.y = element_text(size=12),
						  axis.title=element_text(size=12),
						  legend.position = "none",
						  panel.grid.minor=element_blank(),
						  panel.grid.major.x=element_blank(),
						  plot.margin=unit(c(0,0,0,0.1),"inches"))
			}
		}else{
			X['dID'] <- 'True Bulk'
		}
		if(is.null(totalX)) totalX <- X
		else totalX <- rbind(totalX,X)
	}
	plotX <- plyr::ddply(totalX,c("sID","dID"),function(x)return(data.frame(rmse=sqrt(mean((x$expR-x$estR)^2)))))
	xlabCol <- ifelse(levels(plotX[,"dID"])=="True Bulk","red","grey30")
	xpos <- setNames(1:nlevels(plotX[,'dID']),levels(plotX[,"dID"]))
	xpos["True Bulk"] <- xpos["True Bulk"]+0.7
	plotX <- cbind(plotX,xpos=xpos[as.character(plotX[,"dID"])])
	plotX$sID <- sapply(strsplit(as.character(plotX$sID),"\\|"),function(x)return(paste(x[1:max(2,grep("^GS",x))],collapse="|")))
	
	if(one=="Neuron3"){
		figs[[length(figs)]] <- figs[[length(figs)]]+theme(plot.margin=unit(c(0,0.0,0.4,0.1),"inches"))
		names(xpos) <- gsub("Bulk","Bulk*",names(xpos))
	}
	figs[[length(figs)+1]] <- ggplot(plotX,aes(xpos,rmse,fill=dID))+
		geom_boxplot(outlier.shape=NA,width=bw*nlevels(plotX$dID)/20*(sum(wR)-wR[2])/sum(wR))+
		geom_point(position=position_jitterdodge(),pch=21,size=1.5,color="gray")+
		coord_cartesian(ylim=c(0, ymax[one]))+
		ylab('')+xlab("")+ggtitle(paste(one,"by datasets through each sample"))+
		scale_x_continuous(breaks=xpos,labels=names(xpos))+
		#geom_text_repel(aes(label=sID),data=head(plotX[order(plotX$rmse,decreasing=T),],3),size=2)+
		theme_light()+
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,colour=xlabCol,size=12),
			  axis.text.y = element_text(size=12),
			  axis.title=element_text(size=12),
			  legend.position = "none",
			  panel.grid.minor=element_blank(),
			  panel.grid.major.x=element_blank(),
			  plot.margin=unit(c(0,0.0,0.2,0.1),"inches"))
}
figs <- sapply(figs,function(x)return(list(x+ggtitle(""))))
pdf("Fig4.pdf",width=8.5,height=11)
print(ggarrange(plotlist=figs,labels=LETTERS[1:length(figs)],
				widths=wR,
				ncol=2,nrow=length(figs)/2)+
	  	theme(plot.margin=unit(c(0.5,0.3,0.1,0.3),"inches")))
a <- dev.off()




