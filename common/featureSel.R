###########
## featureSel.R
##
#############
featureSel <- function(X,grp,method="voom",batch=NULL,para,strPrefix=NULL,candiG=NULL){
  saveRDS(X,file="log/cellmapOne.rds")
  cat("Finding features from",nrow(X),"for",length(unique(grp)),"groups by",method,"\n")
  print(table(grp))
  if(!is.null(batch))print(table(batch))
  features <- list()
  for(i in sort(unique(grp))) features[[i]] <- rownames(X)
  strDEG <- NULL
  if(!is.null(strPrefix)) strDEG <- paste(strPrefix,"DESeq.gene.txt",sep=".")
  if(method=="DESeq2"){# for deseq2
    features <- list()
    pheno <- data.frame(row.names=colnames(X),cType=grp)
    if(!is.null(batch)) pheno$batch <- batch
    Data <- DESeqDataSetFromMatrix(countData=matrix(as.integer(as.matrix(X)),nrow=nrow(X),dimnames=dimnames(X)),
                                   colData=pheno,
                                   design=formula(paste("~",paste(colnames(pheno),collapse = "+"))))
    #selG <- apply(X,1,function(x)return(sum(x>16))>2)
    dds <- DESeq(Data,betaPrior=TRUE,parallel=T)#,quiet=T
    oneGrp <- setNames(unique(as.character(grp)),unique(as.character(grp)))
    features<- bplapply(oneGrp,function(i){#
      cat("working on",i,"\n")
      one <- data.frame(row.names=rownames(X))#matrix(0,nrow=nrow(X),ncol=0,dimnames=list(rownames(X),c()))
      for(j in unique(as.character(grp))){
        if(i==j) next
        res <- as.data.frame(results(dds,c("cType",i,j)))
        oneF <- rownames(res)[!is.na(res$padj) & res$padj<para$DEGqvalcut & res$log2FoldChange>para$DEGlogFCcut & res$baseMean>para$DEGbasemeancut]#&selG
        g <- intersect(rownames(one),oneF)
        one <- setNames(cbind(one[g,,drop=F],res[g,"log2FoldChange"]),c(colnames(one),j))
        cat("\tAfter intersect with",j,"DEG:",length(g),"\n")
        #if(length(one)<100) cat("WARNING:",i,"and",j,"are too similar (<100 DEG)\n")
      }
      #print(head(one))
      one <- one[order(apply(one,1,quantile,0.75),decreasing = T),]
      if(nrow(one)<5) cat(paste("Cannot find features for",i))
      return(rownames(one))
    })
    if(F){## one vs rest 
      for(i in unique(grp)){##"Neuron"[sample(9,9)]
        cat("\n\tDESeq2: working on",i,"\n")
        grp1 <- as.character(grp)
        grp1[grp1!=i] <- "others"
        desM1 <- data.frame(row.names=colnames(X),grp=grp1)
        if(!is.null(batch)) desM1 <- cbind(desM1,batch=gsub("\\-","_",batch))
        print(table(grp1))
        Data1 <- DESeqDataSetFromMatrix(countData=matrix(as.integer(as.matrix(X)),nrow=nrow(X),dimnames=dimnames(X)),#counts(dds),#
                                        colData=desM1,
                                        design=formula(paste("~",paste(colnames(desM1),collapse = "+"))))
        dds1 <- DESeq(Data1,betaPrior=TRUE,parallel=T)#,quiet=T
        res <- results(dds1,c("grp",i,"others"))
        rm(list=c("Data1","dds1"))
        if(sum(!is.na(res$padj)&res$padj<0.01&res$log2FoldChange>1)<3){
          cat("\tFew signature genes for",i,"\n")
          next
        }
        res <- res[!is.na(res$padj)&res$padj<0.01,]
        one <- rownames(res)[order(res$log2FoldChange,decreasing=T)[1:sum(res$log2FoldChange>1)]]#min(featureN,)
        cat("\tSignificante UP:",length(one),"\n")
        #cat(length(one),"features selected\n")
        ## remove DEGs among data sets of the same cell type
        cat("\t\tRemove DEG between data sets of the same cell type\n")
        selBatch <- intersect(unique(as.character(batch[grp1!="others"])),c("GSE76381","GSE103723","phs001836","GSE67835","GSE104276",#Brain related cell types
                                                                            "GSE120795","GSE86469","GSE83139"))#pancreatic
        if(length(selBatch)>1){
          for(j in 1:(length(selBatch)-1)){
            for(k in (j+1):length(selBatch)){
              selA <- grp1!="others" & batch==selBatch[j]
              selB <- grp1!="others" & batch==selBatch[k]
              cat("\t\t",selBatch[j],sum(selA),selBatch[k],sum(selB),"\n")
              D1 <- cbind(X[,selA],X[,selB])
              D1 <- D1[apply(D1,1,function(x){return(sum(x>4))})>0.8*min(sum(selA),sum(selB)),]
              compInfo <- data.frame(row.names = colnames(D1),grp=c(rep("A",sum(selA)),rep("B",sum(selB))))
              Data1 <- DESeqDataSetFromMatrix(countData=matrix(as.integer(as.matrix(D1)),nrow=nrow(D1),dimnames=dimnames(D1)),#counts(dds),#
                                              colData=compInfo,
                                              design=formula("~grp"))
              dds1 <- DESeq(Data1,betaPrior=TRUE)#,quiet=T
              res <- results(dds1)
              rm(list=c("D1","Data1","dds1"))
              res <- res[!is.na(res$padj)&res$padj<0.001,]
              one <- one[!one%in%rownames(res)[order(abs(res$log2FoldChange),decreasing = T)][1:min(para$selFeatureN,nrow(res))]]
              cat("\t\tSignificant UP:",length(one),"after compared",selBatch[j],selBatch[k],"\n")
            }
          }
        }
        features[[i]] <- one
      }
    }
    for(i in 1:length(features)){
      uniqueSel <- rep(T,length(features[[i]]))#!features[[i]]%in%unlist(features[-i])
      if(!is.null(candiG)) uniqueSel <- features[[i]]%in%candiG
      features[[i]] <- features[[i]][uniqueSel][1:min(para$selFeatureN,sum(uniqueSel))]
    }
  }else if(method=="Voom"){
    d0 <- DGEList(matrix(as.integer(as.matrix(X)),nrow=nrow(X),dimnames=dimnames(X)))
    d0 <- calcNormFactors(d0)
    d <- d0[apply(cpm(d0),1,function(x){return(sum(x>para$DEGbasemeancut))})>ncol(X)/2,]
    grp <- as.factor(grp)
    mm <- model.matrix(~0 + grp)
    y <- voom(d,mm,plot=T)
    fit <- lmFit(y,mm)
    for(i in levels(grp)){
      ## pair-wised comparisons ----
      one <- c()
      for(j in levels(grp)){
        if(i==j) next
        contr <- makeContrasts(contrasts=paste("grp",i,"-","grp",j,sep=""), levels = colnames(coef(fit)))
        tmp <- eBayes(contrasts.fit(fit, contr))
        res <- topTable(tmp, sort.by = "P", n = Inf,p.value=para$DEGqvalcut)
        if(sum(abs(res$logFC)>1)<100) cat("WARNING:",i,"and",j,"are too similar (<100 DEG)\n")
        one <- c(one,rownames(res)[order(res$logFC,decreasing = T)][1:min(20,sum(res$logFC>para$DEGlogFCcut))])
      }
      if(length(one)<5){
        cat(paste("Cannot find features for",i))
        browser()
      }
      ## one vs all ----------
      grp1 <- as.character(grp)
      grp1[grp!=i] <- "others"
      grp1 <- as.factor(grp1)
      mm1 <- model.matrix(~0 + grp1)
      y1 <- voom(d,mm1,plot=T)
      fit1 <- lmFit(y1,mm1)
      contr1 <- makeContrasts(contrasts=paste("grp1",i,"-","grp1others",sep=""), levels = colnames(coef(fit1)))
      tmp <- eBayes(contrasts.fit(fit1, contr1))
      res <- topTable(tmp, sort.by = "P", n = Inf,p.value=para$DEGqvalcut)
      one <- rownames(res)[order(res$logFC,decreasing = T)][1:min(200,sum(res$logFC>para$DEGlogFCcut))]
      
      features <- unique(c(features,one))
    }
  }else if(method=="Top"){
    features <- list()
    for(i in levels(grp)){
      medianG <- apply(X[,grp==i],1,median)
      features[[i]] <- rownames(X)[order(medianG,decreasing=T)[1:min(para$selFeatureN,sum(medianG>para$DEGbasemeancut))]]
    }
  }else if(method=="edgeR"){
    features <- list()
    dg <- DGEList(2^X-1,group=grp)
    des <- model.matrix(~factor(grp)-1,dg$samples)
    dg <- estimateGLMTrendedDisp(dg,design=des)
    dg <- estimateGLMTagwiseDisp(dg,design=des)
    fit <- glmFit(dg,des)
    gType <- sapply(strsplit(colnames(des),"\\)"),tail,1)
    for(i in gType){
      cat("EdgeR for",i,"\n")
      logFC <- data.frame(row.names=rownames(X))
      for(j in gType){
        if(i==j) next
        one <- setNames(rep(0,length(gType)),gType)
        one[c(i,j)] <- c(1,-1)
        one <- topTags(glmLRT(fit,contrast=one),n=10000)$table
        one <- setNames(one[one$FDR<para$DEGqvalcut&one$logFC>para$DEGlogFCcut,"logFC",drop=F],j)
        logFC <- merge(one,logFC,by="row.names",sort=F)
        rownames(logFC) <- logFC[,1]
        logFC <- logFC[,-1,drop=F]
        cat("\tvs",j,nrow(logFC),"\n")
      }
      features[[i]] <- rownames(logFC)[order(apply(logFC,1,median),decreasing=T)][1:min(nrow(logFC),para$selFeatureN)]
    }
  }else{
    cat("\n\nUNKNOW METHODS: no feature selection\n\n")
  }
  return(features)
}