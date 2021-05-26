require("ComplexHeatmap")
require("circlize")

rm(list=ls())
graphics.off()
closeAllConnections()

Profiles <- c("Major9","CNS6","Neuron3")#,"CNS6","Neuron3"
heatTitle <- setNames(c("log2(TPM+1)","log2(CPM+1)","log2(CPM+1)"),Profiles)
heatCol <- c("#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
			 "#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5",
             "#FFFFFF",
			 "#FFEFEA","#FBE0D1","#FAC9B1","#F5AD8C","#E88B6E","#D96752",
			 "#C6413E","#B31B2C","#8E0C25","#67000d")
#heatCol <- c(#"#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
#        #"#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5",
#        "#FFFFFF",
#        "#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a",
#        "#ef3b2c","#cb181d","#a50f15","#67000d")

pdf("FigS.heatmap.pdf",width=9)
for(one in Profiles){
	D <- readRDS(paste0("../../profiles/",one,".rds"))
	#selG <- unlist(D$selG)
	## order the genes within each cell type profiles
	selG <- sapply(sort(names(D$selG)),function(x){
	    return(D$selG[[x]][order(apply(D$expr[D$selG[[x]],sapply(strsplit(colnames(D$expr),"\\|"),head,1)==x],1,median),
	                        decreasing=T)])
	})
	#browser()
	Exp <- as.matrix(D$expr[unique(unlist(selG)),order(colnames(D$expr))])
	
	## annotation for the gene
	geneDF <- sapply(names(selG),function(x){
		tmp <- rep(x,nrow(Exp))
		tmp[!rownames(Exp)%in%D$selG[[x]]] <- "Other"
		return(tmp)
		})
	geneCol <- sapply(names(D$selG),function(x)return(list(setNames(c(D$para$cellCol[x],"#FFFFFF"),
																	c(x,"Other")))))
	colnames(geneDF) <- names(geneCol) <- paste0(substr(colnames(geneDF),1,3),". genes")
	geneAnno <- HeatmapAnnotation(
		df=data.frame(geneDF,check.names=F),
		col= geneCol,
		#simple_anno_size = unit(0.2, "inch"),
		#annotation_name_side="left",
		show_legend=F
	)
	sampleAnno <- rowAnnotation(
		cellType=sapply(strsplit(colnames(Exp),"\\|"),head,1),
		col=list(cellType=D$para$cellCol),
		show_annotation_name=F,
		show_legend=F
	)
	steps <- quantile(Exp[Exp>0],probs=round((1:(length(heatCol)-0))/(length(heatCol)-0),2))
	#steps <- quantile(Exp[Exp>0],probs=sort(c(round((1:(length(heatCol)-2))/(length(heatCol)-2),2),0.95,0.98)))
	#print(steps)
	#steps <- seq(0,max(Exp),length.out=length(heatCol))
	col_fun <- circlize::colorRamp2(steps,heatCol)
	
	ht <-Heatmap(t(Exp),name=heatTitle[one],col=col_fun,
			show_row_names=F,show_column_names=F,
			cluster_rows=F,cluster_columns=F,
			column_title="Profile genes",#row_title="Profile cell type",
			row_split = paste0(substr(sapply(strsplit(colnames(Exp),"\\|"),head,1),1,3),". profiles"),
			row_title_rot=0,row_title_gp=gpar(fontsize=12),
			top_annotation=geneAnno,left_annotation=sampleAnno,
			heatmap_legend_param=list(legend_height = unit(2, "inches"),
									  title_position = "lefttop-rot"),
			use_raster = TRUE)
	draw(ht)
}
a <- dev.off()