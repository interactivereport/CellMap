require(ggplot2)
require(ggpubr)
require(reshape2)
require(tidyverse)

figData = readRDS('Fig6.data.rds')

#functions
calculate_RMSE <- function(decomp, ground){
    ground <- ground[,colnames(decomp)]
    sqrd_diff <- (ground - decomp)^2
    RMSE <- sqrt(colSums(sqrd_diff)/nrow(sqrd_diff))
    return(RMSE)
}

plot_scatters_mixed <- function(tpm, ground, title){

    sorted_cols <- sort(colnames(ground))
    ground<-ground[,sorted_cols]
    tpm<-as.data.frame(tpm[,sorted_cols])
    ground_melt <- melt(as.matrix(ground))
    colnames(ground_melt) <- c('Cell_Type','Sample','Ground')
    edgeR_tpm_melt <- melt(tpm)
    colnames(edgeR_tpm_melt) <- c('Sample2','CellMap')
    
    rmse <- mean(calculate_RMSE(tpm,ground))

    scatter_data <- cbind(ground_melt,edgeR_tpm_melt)
    scatter_data$Ground <- scatter_data$Ground 
    scatter_data$CellMap <- scatter_data$CellMap 
    scatter_data$difference <- scatter_data$CellMap - scatter_data$Ground 

    avg_data <- scatter_data %>% select(Cell_Type,Ground,CellMap) %>% group_by(Cell_Type) %>% summarize(Ground_SD = sd(Ground),CellMap_SD=sd(CellMap),Ground=mean(Ground),CellMap=mean(CellMap))
    
    p <- ggplot(scatter_data, aes(y=CellMap, x=Ground, color=Cell_Type)) + 
        xlim(0,1)+
        ylim(0,1)+
        xlab(paste(title,'Ground Truth'))+
        ylab('CellMap Prediction')+
        geom_point(size=3, shape=1, alpha=.4)+
        geom_abline(intercept = 0, slope = 1, color="black", 
                     linetype="dashed", size=.5)+
        geom_point(data=avg_data, size=4.5, aes(color=Cell_Type))+
        geom_errorbarh(data=avg_data, aes(y=CellMap, xmin=pmax(0,Ground - Ground_SD), xmax=Ground + Ground_SD), alpha=.9,height =.02)+
        geom_errorbar(data=avg_data, aes(x=Ground, ymin=pmax(0,CellMap - CellMap_SD), ymax=CellMap + CellMap_SD), alpha=.9,width=.02)+ 
        theme_light()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              aspect.ratio = 1,
              legend.position = "none",
              panel.grid.minor=element_blank())
    return(p)
}

# A cadio -------
adf <- figData$A

cadio <- ggplot(adf,aes(x=time,y=Cardiomyocytes,fill=Cardiomyocytes))+
	#geom_boxplot(width=.5,fill='#99445E',position = position_dodge2(preserve = "single"))+
    geom_boxplot(outlier.shape=NA,fill='#99445E')+
    geom_point(position=position_jitterdodge(),size=1)+
	ylab('Est. Cardio. (%)')+xlab("")+
	#geom_jitter(alpha=0.7,size=0.5,width=.05)+
	theme_light()+
	#scale_y_continuous(breaks=c(98,99,100))+
	theme(axis.text.x=element_text(angle=45,hjust=0.8,vjust=0.9,size=12),
		  axis.title=element_text(size=11),
		  panel.grid.minor=element_blank(),
		  panel.grid.major.x=element_blank(),
		  legend.position = "none",
		  plot.margin=unit(c(0.3,0.1,0,0),"inches"))

# B neuron--------
B <- figData$B

neuronC <- ggplot(B$estC,aes(x=time,y=neuron,fill=batch))+
	geom_boxplot(color='gray')+
	ylab('Est. neuron (%)')+labs(fill="Line")+xlab("")+
	geom_jitter(alpha=0.7,size=1,position=position_jitterdodge(jitter.width=0.2))+
	theme_light()+
	scale_y_continuous(breaks=c(92,94,96,98,100))+
	theme(axis.text.x=element_text(size=12),
		  axis.title=element_text(size=11),
		  panel.grid.minor=element_blank(),
		  panel.grid.major.x=element_blank(),
		  legend.position = "bottom",
		  plot.margin=unit(c(0.3,0.1,0,0),"inches"))

neuronP <- ggplot(B$estP,aes(x=time,y=neuron,color=batch,shape=batch))+
	geom_point()+
	stat_smooth(aes(group=batch),se=F,method ='loess',formula=y~x,size=1)+
	ylab('-log10(p-value)')+xlab("Time")+
	theme_light()+
	theme(axis.text.x=element_text(size=12),
		  axis.title=element_text(size=11),
		  panel.grid.minor=element_blank(),
		  panel.grid.major.x=element_blank(),
		  legend.position="none",
		  plot.margin=unit(c(0.3,0.1,0,0),"inches"))

#neuron <- ggarrange(fD+rremove('x.text')+rremove("xlab"),fP,ncol=1,nrow=3,heights = c(1,1,1))+
#    theme(plot.margin=unit(c(0.1,0.1,0,0),"inches"))

# C Microglia---------
C <- figData$C

micro <- ggplot(C$estC,aes(x=time,y=microglia,fill=treatment))+
	geom_boxplot(aes(color=treatment),width=1,position = position_dodge2(preserve = "single"))+
	ylab('Est. microglia (%)')+xlab("")+
    geom_jitter(alpha=0.4,size=1,position=position_dodge(1))+
	theme_light()+
	#scale_y_continuous(breaks=c(98,99,100))+
	theme(axis.text.x=element_text(angle=45,hjust=0.8,vjust=0.9,size=12),
		  axis.title=element_text(size=12),
		  panel.grid.minor=element_blank(),
		  panel.grid.major.x=element_blank(),
		  legend.position = c(0.8, 0.4),
		  legend.text = element_text(size=10),
		  legend.key.size=unit(0.8,'lines'),
		  plot.title = element_text(size = 13),
		  plot.margin=unit(c(0.3,0.1,0,0),"inches"))+
    guides(fill=guide_legend(override.aes = list(size=1)))

# D ROSMAP -----------
D = figData$D
cols <- readRDS(paste0("../../profiles/Major9.rds"))$para$cellCol
ros_decomp = D$ros_decomp
rosIHC_decomp = D$rosIHC_decomp
rossnRNA_decomp = D$rossnRNA_decomp

# IHC vs CM
for(row in rownames(ros_decomp)){
    if(!row %in% rownames(rosIHC_decomp)){
        rosIHC_decomp[row,] <- 0
    }
}
rosIHC_decomp <- rosIHC_decomp[order(rownames(rosIHC_decomp)),]
rosIHC_decomp <- t(t(rosIHC_decomp) / colSums(rosIHC_decomp))
selected_cols <- colnames(ros_decomp)[colnames(ros_decomp) %in% colnames(rosIHC_decomp)]
ros_decomp_sub <- ros_decomp[selected_cols]
ros_decomp_sub <- ros_decomp_sub[order(colnames(ros_decomp_sub))]
rosIHC <- plot_scatters_mixed(ros_decomp_sub,rosIHC_decomp,'IHC') + 
    coord_flip()+
    scale_color_manual(values=cols)+
    theme(legend.position="bottom",
          legend.title=element_blank(),
          legend.text = element_text(size=9),
          legend.key.size=unit(0.3,'lines'),
          plot.margin=unit(c(0.25,0,0.1,0.1),"inches"))+
    guides(color=guide_legend(ncol=3,override.aes = list(size=2)))

# snRNA vs CM
missing_cols <- rownames(ros_decomp)[!rownames(ros_decomp) %in% rownames(rossnRNA_decomp)]
rossnRNA_decomp[missing_cols,] <- 0
selected_cols <- colnames(rossnRNA_decomp)[colnames(rossnRNA_decomp) %in% colnames(ros_decomp)]
rossnRNA_decomp <- rossnRNA_decomp[order(rownames(rossnRNA_decomp)),selected_cols]
selected_cols <- colnames(ros_decomp)[colnames(ros_decomp) %in% colnames(rossnRNA_decomp)]
ros_decomp_sub <- ros_decomp[selected_cols]
ros_decomp_sub <- ros_decomp_sub[order(colnames(ros_decomp_sub))]
rosRNA <- plot_scatters_mixed(ros_decomp_sub,rossnRNA_decomp,'snRNA') +
    coord_flip()+
    scale_color_manual(values=cols)+
    theme(legend.position="none",
          plot.margin=unit(c(0.25,0,0.1,0.1),"inches"))

# # IHC vs snRNA
# ground1 <- rosIHC_decomp 
# ground2 <- rossnRNA_decomp
# select <- colnames(ground1) %in% colnames(ground2)
# select2 <- colnames(ground2) %in% colnames(ground1)
# ground1 <-ground1[rownames(ground2),select]
# ground2 <-ground2[,select2]
# Dp3 <- plot_scatters_mixed(ground1,ground2,'ROSMAP snRNA') + ylab('IHC Ground Truth') + 
#     xlab('snRNA Ground Truth') + theme(legend.position = "none") + scale_color_manual(values=cols)
#ros <- ggarrange(Dp1+rremove('x.text')+rremove("xlab"),Dp2,ncol=1,nrow=2,heights = c(1.05,1))+
#    theme(plot.margin=unit(c(0,0,0,0),"inches"))

# saving cadio, neuronC, neuronP, micro, rosIHC, rosRNA---------
left <- ggarrange(cadio,micro,neuronC,neuronP,ncol=1,nrow=4,labels = c("A", "B", "C", "D"),heights = c(1,1.1,1.2,1))
right <- ggarrange(rosIHC,rosRNA,ncol=1,nrow=3,labels = c("E", "F"),heights = c(1.05,0.9,0.2))


pdf('Fig6.pdf',width=8.5,height=11)
print(ggarrange(left,right,widths = c(1,1))+
          theme(plot.margin=unit(c(0.5,0.3,0.3,0.3),"inches")))
a <- dev.off()
