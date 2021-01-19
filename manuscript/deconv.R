##############################################
## In order for the following to run, the compiled all single cells expression of all correspond public data sets is required.
## An RDS file contains expression matrix for each of Major9, CNS6, Neuron3 profile should located in defined 'strProfile' folder.
## The rows of the matrix are the genes (row names), while the columns are cells (each of the column names consist of cellType|dataset|...).
## Please make sure the cell names (column names) are the correct format (above) to include cell type and dataset information.
##
#############################################
rm(list=ls())
gc()
graphics.off()
closeAllConnections()
strRoot <- "./"
cellmapProfile <- paste0(strRoot,"profiles/")
strProfile <- paste0(strRoot,"manuscript/profiles/")
strPath <- paste0(strRoot,"manuscript/pseudo/")
strOut <- paste0(strRoot,"manuscript/deconv/")

###
loadDECONV <- function(){
	source(paste0(strRoot,"manuscript/common/deconv.Bisque.R"))
	source(paste0(strRoot,"manuscript/common/deconv.SCDC.R"))
	source(paste0(strRoot,"manuscript/common/deconv.MuSiC.R"))
	source(paste0(strRoot,"manuscript/common/deconv.CellMap.R"))
	source(paste0(strRoot,"manuscript/common/plotEval.R"))
}
suppressWarnings(suppressMessages(loadDECONV()))
methods <- c("CellMap","MuSiC","Bisque","SCDC")#,
strD <- c("Major9","CNS6","Neuron3")#

if(T){
	register(MulticoreParam(32))
	tRMSE <- list()
	for(dID in strD){
		for(dType in c("pure","mix")){
			message("===============",dID," ",dType,"================")
			bulk <- read.table(paste0(strPath,dID,".",dType,".txt"),header=T,row.names=1,sep="\t",as.is=T,check.names=F)
			expR <- read.table(paste0(strPath,dID,".",dType,".rate"),header=T,row.names=1,sep="\t",as.is=T,check.names=F)
			expR <- apply(expR,2,function(x)return(x/sum(x)))
			rmse <- c()
			for(m in methods){
				message("---------------------",m,"---------------------")
				strF <- paste0(strOut,dID,".",dType,".",m,".rdata")
				if(file.exists(strF)){
					load(strF)
				}else{
					if(m=="CellMap"){
						Profi <- readRDS(paste0(cellmapProfile,dID,".rds"))
					}else{
						Profi <- readRDS(paste0(strProfile,dID,".rds"))
					}
					#save(bulk,Profi,file=gsub("rdata$","tmp.rdata",strF))
					print(system.time({cellComp <- get(paste0("deconv.",m))(bulk,Profi)}))
					cellCol <- Profi$para$cellCol
					save(cellComp,cellCol,file=strF)
				}
				estR <- cellComp$composition
				estR[cellComp$compoP>0.05] <- 0
				estR <- estR[rownames(expR),]
				pdf(gsub("rdata$","pdf",strF))
				rmse <- cbind(rmse,plotEval(estR,expR,cellCol)[,"RMSE"])
				dev.off()
			}
			dimnames(rmse) <- list(colnames(bulk),methods)
			tRMSE[[paste(dID,dType,sep="|")]] <- rmse
			saveRDS(tRMSE,file=paste0(strOut,"rmse.rds"))
		}
	}
}
tRMSE <- readRDS(paste0(strOut,"rmse.rds"))
yMax <- setNames(c(0.5,0.5,0.6,0.6,0.8,0.8),names(tRMSE))
pdf(paste0(strOut,"overall.pdf"))
for(one in names(tRMSE)){
	D <- melt(tRMSE[[one]])
	colnames(D) <- c("sID","method","RMSE")
	X <- ddply(D,'method',summarise,mu=mean(RMSE),sigma=sd(RMSE))
	p <- ggplot(X,aes(method,mu,fill=method))+
		geom_bar(stat='identity',color='black')+
		geom_point(aes(method,RMSE),data=D,position=position_jitterdodge(),pch=21,size=1,color="gray")+
		geom_errorbar(aes(ymin=mu,ymax=mu+sigma),width=0.3)+
		ylim(0,yMax[one])+
		#scale_fill_manual(values=col)+
		ylab('RMSE')+xlab("")+ggtitle(one)+
		theme_classic()+
		theme(axis.text.x = element_text(angle=90,hjust=1))
	print(p)
	
}
dev.off()