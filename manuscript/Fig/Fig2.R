require(pheatmap)
require(RColorBrewer)

rm(list=ls())
graphics.off()
closeAllConnections()

Dinfo <- readRDS("Fig2.rds")
heatCol <- c("#FFFFFFFF",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(19))


pdf('Fig2.pdf',width=15)
for(one in names(Dinfo)){
	message("===========",one,"================")
	D <- Dinfo[[one]]
	para <- readRDS(paste0("../../profiles/",one,".rds"))$para
	
	pheatmap(D,heatCol,breaks=2^seq(0,ceiling(log2(max(D))),length.out=21),
			 cellwidth=40,cellheight=30,number_format = "%d",
			 main=one,
			 cluster_rows=F,cluster_col=F,legend=F,display_numbers=T,fontsize_number=10,
			 annotation_col=para$tech[colnames(D),,drop=F])
}
dev.off()








