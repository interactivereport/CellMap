require(ggpubr)
require(ggplot2)
require(reshape2)

rm(list=ls())
graphics.off()
closeAllConnections()

tDiff <- readRDS("S3.4.5.rds")
methods <- c("CellMap",'MuSiC','Bisque','SCDC')
avgData <- F

Profiles <- c(Major9=3,CNS6=4,Neuron3=5)#
ymax <- setNames(c(0.5,0.7,0.7),names(Profiles))
bw <- 2.2
wR <- c(1,1,1,1)

for(one in names(Profiles)){
	figs <- list()
	cellCol <- readRDS(paste0("../../profiles/",one,".rds"))$para$cellCol
	for(pseudo in c("pure","mix")){
		X <- tDiff[[paste(one,pseudo,sep='|')]]
		for(m in methods){
			if(avgData){
				plotX <- plyr::ddply(X,c("celltype"),function(x)return(data.frame(rmse=sqrt(mean((x$expR-x[,m])^2)))))
				figs[[length(figs)+1]] <- ggplot(plotX,aes(celltype,rmse,fill=celltype))+
					geom_bar(stat="identity")+
					coord_cartesian(ylim=c(0, ymax[one]))
			}
			else{
				X['dID'] <- sapply(strsplit(as.character(X$sID),"\\|"),"[[",2)
				plotX <- plyr::ddply(X,c("celltype","dID"),function(x)return(data.frame(rmse=sqrt(mean((x$expR-x[,m])^2)))))
				#labelX <- plotX[NULL,]
				#for(i in levels(plotX$celltype))labelX <- rbind(labelX,head(plotX[intersect(order(plotX$rmse,decreasing=T),grep(i,plotX$celltype)),],1))
				#labelX <- labelX[labelX$rmse>0.05,]
				#plotX$celltype <- factor(plotX$celltype,levels=names(cellCol))
				figs[[length(figs)+1]] <- ggplot(plotX,aes(celltype,rmse,fill=celltype))+
					geom_boxplot(outlier.shape=NA,width=0.8*length(cellCol)/9)+#bw*nlevels(plotX$celltype)/20*(sum(wR)-wR[1])/sum(wR)
					geom_point(position=position_jitterdodge(),pch=21,size=1.5,color="gray")+
					#geom_jitter(alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.2))+
					coord_cartesian(ylim=c(0, ymax[one]))
					#geom_text_repel(aes(label=dID),data=labelX,size=2)+
			}
			figs[[length(figs)]] <- figs[[length(figs)]]+
				scale_fill_manual(values=cellCol)+
				ylab("")+xlab("")+ggtitle(paste(m,pseudo,sep=" "))+
				theme_light()+
				theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=12),
					  axis.text.y = element_text(size=12),
					  axis.title=element_text(size=12),
					  legend.position = "none",
					  panel.grid.minor=element_blank(),
					  panel.grid.major.x=element_blank(),
					  plot.margin=unit(c(0.2,0,0,0),"inches"))
			if(m=='CellMap') figs[[length(figs)]] <- figs[[length(figs)]] + ylab("RMSE")
			
		}
	}
	if(avgData)pdf(paste0("S",Profiles[one],"_.pdf"),width=8.5,height=11)
	else pdf(paste0("S",Profiles[one],".pdf"),width=8.5,height=11)
	print(ggarrange(plotlist=figs,labels=LETTERS[1:length(figs)],
					ncol=4,nrow=3)+
		  	theme(plot.margin=unit(c(0.3,0.2,0.1,0.2),"inches")))
	a <- dev.off()
}



