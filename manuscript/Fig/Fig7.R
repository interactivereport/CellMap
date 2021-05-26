#require(plyr)
require(ggpubr)
require(ggplot2)
require(reshape2)

rm(list=ls())
graphics.off()
closeAllConnections()


tRMSE <- readRDS("Fig7.rds")
tRMSE <- tRMSE[c(grep('pure$',names(tRMSE)),grep('mix$',names(tRMSE)))]
yMax <- setNames(c(0.8,0.8,0.8,0.8,0.8,0.8),names(tRMSE))
figs <- list()
for(one in names(tRMSE)){
	D <- melt(tRMSE[[one]])
	colnames(D) <- c("sID","method","RMSE")
	X <- plyr::ddply(D,'method',summarise,mu=mean(RMSE),sigma=sd(RMSE))
	figs[[length(figs)+1]] <- ggplot(X,aes(method,mu,fill=method))+
		geom_bar(stat='identity',color='black')+
		geom_point(aes(method,RMSE),data=D,position=position_jitterdodge(),pch=21,size=1,color="gray")+
		geom_errorbar(aes(ymin=mu,ymax=mu+sigma),width=0.3)+
		ylim(0,yMax[one])+
		xlab("")+ylab(paste0('RMSE (',one,')'))+#+ggtitle(paste("Profile:",one))
		theme_light()+
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=12),
			  axis.text.y = element_text(size=12),
			  axis.title=element_text(size=12),
			  legend.position = "none",
			  panel.grid.minor=element_blank(),
			  panel.grid.major.x=element_blank(),
			  plot.margin=unit(c(0,0,0,0.1),"inches"))
	
}
figs <- sapply(figs,function(x)return(list(x+ggtitle(""))))
pdf("Fig7.pdf",width=8.5,height=11)
print(ggarrange(plotlist=figs,labels=LETTERS[1:length(figs)],
				ncol=3,nrow=length(figs)/2)+
	  	theme(plot.margin=unit(c(0.5,0.3,0.1,0.3),"inches")))
a<- dev.off()