require(reshape2)
require(ggpubr)
require(ggplot2)
require(ggrepel)
require(tidyverse)

#functions
create_gt_df <- function(df, ctypes){
    cnms <- colnames(df)
    rnms <- rownames(df)
    n <- length(rnms)
    for(i in c(1:length(cnms))){
        col <- cnms[i]
        ct <- ctypes[i]
        
        if(ct=='Skip'){
            values = rep(NA,n)
        } else{
            values = rep(0,n)
        }
        
        if(i==1){
            df <- data.frame(a = values)
            colnames(df) <- col
            rownames(df) <- rnms
        } else {
            df[col] = values
        }
        
        if(ct!='Skip'){
            df[ct,col] <- 1
        }
    }
    return(df)
}

calculate_RMSE <- function(decomp, ground){
    ground <- ground[,colnames(decomp)]
    sqrd_diff <- (ground - decomp)^2
    RMSE <- sqrt(colSums(sqrd_diff)/nrow(sqrd_diff))
    return(RMSE)
}


evals_list <- readRDS('Fig5.data.rds')

#evals_list <- list(CNS6=CNS6_evals,Major9=Major9_evals) -----
colors <- c('#e41a1c', # 
    '#377eb8', #
    '#984ea3', #
    '#a65628', #
    '#ff7f00', #
    '#2a1286') #

# ------
plots <- list()
i <- 1
for(eval_name in rev(names(evals_list))){
    colors <- readRDS(paste0("../../profiles/",eval_name,".rds"))$para$cellCol
    evals <- evals_list[[eval_name]]
    df_ct <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(df_ct) <- c('celltype','RMSE','reference')

    df_data <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(df_data) <- c('Dataset','RMSE')

    for(e_eval in evals){
        type <- e_eval$type
        reference <- gsub('Pure','Primary',e_eval$reference)
        if(type %in% c('ROSMAP-snRNA','ROSMAP-IHC')){
            next
        }

        # create data table for RMSE by CT
        RMSE <- e_eval$RMSE
        df_ct <- rbind(df_ct,data.frame(celltype=type,RMSE=RMSE,reference=reference))

        # create data table for RMSE by dataset
        datasets <- e_eval$data
        dnames <- names(datasets)
        for(name in dnames){
            ds <- datasets[[name]]
            ds_gt <- create_gt_df(ds,rep(type,dim(ds)[2])) #calculate the ground truth table
            ds_RMSE <- calculate_RMSE(ds,ds_gt)
            ds_name <- rep(name,length(ds_RMSE))
            ds_type <- rep(type,length(ds_RMSE))
            ds_reference <- rep(reference,length(ds_RMSE))
            
            df_data <- rbind(df_data,data.frame(Dataset=ds_name,RMSE=ds_RMSE,celltype=ds_type,reference=ds_reference))
        }
    }
    
    df_median <- df_data %>% group_by(Dataset,celltype,reference) %>% summarize(RMSE=median(RMSE))
    
    df_ct$reference <- gsub('Pure','Purified_Primary',df_ct$reference)
    df_ct$reference <- factor(df_ct$reference,levels = c('Purified_Primary','iPSC'))
    
    bCol <- 0.1
    fillColors <- colors[levels(df_data$celltype)]
    fillAlpha <- setNames(c(1,0),levels(df_data$reference))
    df_data$cellRef <- paste(df_data$celltype,df_data$reference,sep="|")
    borderColor <- setNames(rep(fillColors,length(fillAlpha)),
                            paste(rep(names(fillColors),length(fillAlpha)),
                                  sort(rep(names(fillAlpha),length(fillColors))),sep="|"))
    borderColor[grepl(names(fillAlpha)[fillAlpha==1],names(borderColor))] <- rgb(bCol,bCol,bCol)
    legendFill <- setNames(rgb(bCol+0.4,bCol+0.4,bCol+0.4,alpha=fillAlpha),
                           names(fillAlpha))
    legendCol <- setNames(rep(rgb(bCol+0.4,bCol+0.4,bCol+0.4),length(fillAlpha)),names(fillAlpha))
    legendCol[fillAlpha==1] <- rgb(bCol,bCol,bCol)
    
    
    
    plots[[i]]<-ggplot(df_data,aes(celltype,RMSE,fill=celltype))+ #,color=cellRef
        facet_wrap(~reference,ncol=2)+#
        geom_boxplot(outlier.shape=NA)+ #,aes(alpha=reference)
        geom_jitter(width=0.2,pch=21,size=1.5,color='gray')+ #aes(group=reference),
        ylab('RMSE')+xlab("")+labs(alpha="Origin:")+
        scale_fill_manual(values=fillColors)+
        coord_cartesian(ylim=c(0, 0.5))+
        scale_color_manual(values=borderColor)+
        scale_alpha_manual(values=fillAlpha)+
        theme_light()+
        theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=10),
              axis.text.y = element_text(size=12),
              axis.title=element_text(size=12),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              strip.text = element_text(size = 10,color="black"),
              strip.background = element_blank(), #element_rect(fill="white",colour="gray"),
              plot.margin=unit(c(0.2,0,0.1,0.1),"inches"))+
        guides(fill=F,color=F,
               alpha=guide_legend(override.aes = list(fill=legendFill,color=legendCol)))

    #dataset_order <- unique(c(sort(unique((as.character(df_data[df_data$reference==levels(df_data$reference)[1],]$Dataset)))) ,
    #                          sort(unique((as.character(df_data[df_data$reference==levels(df_data$reference)[2],]$Dataset))))))
    #df_data$Dataset <- factor(as.character(df_data$Dataset),levels = dataset_order)
    #legendFill <- paste0("#4daf4a",substr(legendFill,8,9))
    #legendCol[2] <- "#4daf4a"
    
    ## change the in-house data:
    #dL <- levels(df_data$Dataset)
    #dL[grepl("GSEXXXX",dL)] <- paste0(gsub("GSEXXXX_","",grep("GSEXXXX",dL,value=T)),"*")
    #ds <- as.character(df_data$Dataset)
    #ds[grepl("GSEXXXX",ds)] <- paste0(gsub("GSEXXXX_","",grep("GSEXXXX",ds,value=T)),"*")
    #df_data$Dataset <- factor(ds,levels=dL)
    ds <- gsub('_adult|GSEXXXX_|_iPSC','',as.character(df_data$Dataset))
    ds[grepl('Batch',ds)] <- paste0(ds[grepl('Batch',ds)],"*")
    df_data$Dataset <- ds
    
    plots[[i+1]]<-ggplot(df_data,aes(x=Dataset,y=RMSE))+#,fill=reference,color=reference
                        geom_boxplot(outlier.shape=NA,fill='#4daf4a')+
                        geom_jitter(width=0.2,pch=21,size=1.5,color='gray',fill='black')+
                        #geom_jitter(alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.2))+
                        facet_grid(.~reference,scales = 'free',space='free')+#
                        coord_cartesian(ylim=c(0, 0.5))+
                        ylab('RMSE')+xlab("")+
                        scale_fill_manual(values=legendFill)+
                        scale_color_manual(values=legendCol)+
                        theme_light()+
                        theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=10),
                              axis.text.y = element_text(size=12),
                              axis.title=element_text(size=12),
                              legend.position = "none",
                              panel.grid.minor=element_blank(),
                              panel.grid.major.x=element_blank(),
                              strip.text = element_text(size = 10,color="black"),
                              strip.background = element_blank(), #element_rect(fill="white",colour="gray"),
                              plot.margin=unit(c(0.1,0,0,0.1),"inches"))
    i <- i + 2
}

## arrange the plots

AB <- ggarrange(plotlist=list(plots[[1]],plots[[3]]),ncol=2,nrow=1,labels = c("A", "B"),
                common.legend = TRUE,legend="bottom")+
    theme(plot.margin=unit(c(0,0,0,0),"inches"))


# sort(unique((as.character(df_data[df_data$reference=='Pure',]$Dataset)))))
pdf('Fig5.pdf',width=8.5,height=11)
print(ggarrange(plotlist=list(AB,plots[[2]],plots[[4]]),ncol=1,nrow=3,labels = c("", "C", "D"))+
          theme(plot.margin=unit(c(0.5,0.3,0.1,0.3),"inches")))
a <- dev.off()