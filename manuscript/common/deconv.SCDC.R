##########################
## deconv.SCDC.R
##
##########################
loadSCDC <- function(){
  if(!require(Biobase)){
    BiocManager::install("Biobase")
    if(!require(Biobase)) stop("Cannot install Biobase!")
  }
  if(!require(L1pack)){
    install.packages("L1pack",repos="https://cloud.r-project.org/")
    if(!require(L1pack)) stop("Cannot install L1pack!")
  }
  if(!require(nnls)){
    install.packages("nnls",repos="https://cloud.r-project.org/")
    if(!require(nnls)) stop("Cannot install nnls!")
  }
  if(!require(lsei)){
    install.packages("lsei",repos="https://cloud.r-project.org/")
    if(!require(lsei)) stop("Cannot install lsei!")
  }
  if(!require(BiocParallel)){
    install.packages("BiocParallel",repos="https://cloud.r-project.org/")
    if(!require(BiocParallel)) stop("Cannot install BiocParallel!")
  }
}
suppressWarnings(suppressMessages(loadSCDC()))

deconv.SCDC <- function(bulk,feature,selCellType=NULL){
  SC <- feature$SC
  rm(feature)
  if(!is.null(selCellType)) SC <- SC[,sapply(strsplit(sampleNames(SC),"\\|"),head,1)%in%selCellType]
  cellSet <- sapply(strsplit(sampleNames(SC),"\\|"),function(x){return(paste(x[1:2],collapse="|"))})
  uCellSet <- unique(cellSet)
  uPairs <- data.frame(cType=sapply(strsplit(uCellSet,"\\|"),head,1),
  					 dID=sapply(strsplit(uCellSet,"\\|"),tail,1),
  					 stringsAsFactors=F)
  cType <- unique(uPairs[,'cType'])
  dID <- unique(uPairs[,'dID'])
  dSets <- list()
  #dProb <- table(sapply(strsplit(cellSet,"\\|"),tail,1))
  dProb <- setNames(rep(10,length(dID)),dID)#(sum(dProb)-dProb)/sum(dProb)#
  while(length(dSets)<3||sum(!dID%in%unlist(dSets))>0){
  	one <- sort(sample(dID,sample(ceiling(length(dID)/1.5),1),prob=dProb))
  	#cat("Select:",one)
  	selType <- table(uPairs[uPairs[,'dID']%in%one,'cType'])
  	if(length(selType)<length(cType)||min(selType)<2){
  		#cat("\tfailed\n")
  		next
  	}
  	bEx <- F
  	for(a in dSets){
  		if(length(a)==length(one) && sum(a!=one)==0){
  			bEx <- T
  			break
  		}
  	}
  	if(bEx){
  		#cat("\tExisted\n")
  		next  		
  	}
  	cat("Select:",one)
  	cat("\tsuccessful\n")
  	dSets[[length(dSets)+1]] <- one
  	#dProb[one] <- sapply(one,function(x)return(max(1,dProb[x]-1)))
  	dProb[one] <- dProb[one]/2
  }
  message("finished data sets combination with ",length(dSets)," sets")
  sc.eset.list <<- list()
  for(one in dSets){
  	sc.eset.list[[length(sc.eset.list)+1]] <<- SC[,sapply(strsplit(cellSet,"\\|"),tail,1)%in%one]
  }
  
  rm(SC,cellSet)
  gc()
  res <- SCDC_ENSEMBLE(bulk.eset=ExpressionSet(as.matrix(bulk)),ct.varname="cellType",sample="dID",ct.sub=cType)
  cat("\n\n")
  comp <- t(wt_prop(res$w_table[1,1:length(sc.eset.list)],res$prop.only))
  pV <- wt_extP(res$w_table[1,1:length(sc.eset.list)],res$prop.list)
  return(list(composition=comp,
              compoP=head(pV,-1),
              overallP=tail(pV,1),
              coverR=wt_cover(res$prop.list,rownames(sc.eset.list[[1]])),
              rawComp=res$prop.only,
              rawSets=NULL,
              missingF=""))
}

SCDC_ENSEMBLE <- function (bulk.eset,  ct.varname, sample, #sc.eset.list = NULL,
                           ct.sub, grid.search = F, search.length = 0.05, iter.max = 2000, 
                           nu = 1e-04, epsilon = 0.001, truep = NULL, weight.basis = T, 
                           prop.input = NULL, Transform_bisque = F, ...) {
  if (!is.null(prop.input)) {
    message("Using user-input estimated proportion list ...")
    prop.list <- prop.input
  }
  else {
  	totalN <- length(sc.eset.list)
  	iCount <- 1
    prop.list <- lapply(sc.eset.list, function(zz) {#bp
      message(iCount,"/",totalN)
      iCount <<- iCount+1
      if (length(unique(zz@phenoData@data[, sample])) > 1) {
        SCDC_prop(bulk.eset = bulk.eset, sc.eset = zz, 
                  ct.varname = ct.varname, sample = sample, truep = truep, 
                  ct.sub = ct.sub, iter.max = iter.max, nu = nu, 
                  epsilon = epsilon, weight.basis = weight.basis, 
                  Transform_bisque = Transform_bisque)
      }
      else {
        SCDC_prop_ONE(bulk.eset = bulk.eset, sc.eset = zz, 
                      ct.varname = ct.varname, sample = sample, truep = truep, 
                      ct.sub = ct.sub, iter.max = iter.max, nu = nu, 
                      epsilon = epsilon, weight.basis = weight.basis)
      }
    })
  }
  row.list <- sapply(1:length(prop.list), function(x) {
    rownames(prop.list[[x]]$yhat)
  })
  gene.prop <- Reduce("intersect", row.list)
  gene.prop2 <- intersect(gene.prop, rownames(bulk.eset))
  subj.order <- colnames(bulk.eset)
  ycpm <- getCPM0(exprs(bulk.eset)[gene.prop2, subj.order])
  g.filter <- rowSums(ycpm) < quantile(rowSums(ycpm), 0.95) & 
    rowSums(ycpm) > quantile(rowSums(ycpm), 0.15)
  gene.use <- gene.prop2[g.filter]
  length(gene.use)
  yv <- c(getCPM0(exprs(bulk.eset)[gene.use, subj.order])) * 
    1e+05
  y.list <- do.call(cbind, lapply(prop.list, function(x) {
    c(getCPM0(x$yhat[gene.use, subj.order])) * 1e+05
  }))
  sse <- function(x, y) {
    sum((x - y)^2, na.rm = T)
  }
  sae <- function(x, y) {
    sum(abs(x - y), na.rm = T)
  }
  rmsd <- function(x, y) {
    sqrt(mean((x - y)^2, na.rm = T))
  }
  message("Searching ENSEMBLE weight by Sum of Squared Errors or Sum of Abs Errors ......")
  sses <- apply(y.list, 2, function(x) {
    sse(yv, x)
  })
  sse.wt <- 1/sses/sum(1/sses)
  saes <- apply(y.list, 2, function(x) {
    sae(yv, x)
  })
  sae.wt <- 1/saes/sum(1/saes)
  rmsds <- apply(y.list, 2, function(x) {
    rmsd(yv, x)
  })
  rmsd.wt <- 1/rmsds/sum(1/rmsds)
  message("Searching ENSEMBLE weight by LAD -- Minimizing mAD of Y measurement")
  w_lad <- NA
  dt <- data.frame(y = yv, y.list)
  fitlad <- L1pack::lad(y ~ . - 1, data = dt, method = c("BR", 
                                                         "EM"))
  w_lad <- fitlad$coefficients
  w_lad[w_lad < 0] <- 0
  w_lad <- w_lad/sum(w_lad)
  w_nnls <- NA
  message("Searching ENSEMBLE weight by NNLS -- Minimizing MSE of Y measurement")
  fitnnls <- nnls::nnls(A = as.matrix(y.list), b = yv)
  w_nnls <- fitnnls$x
  w_nnls <- w_nnls/sum(w_nnls)
  combo <- lapply(prop.list, function(x) {
    x$prop.est.mvw
  })
  w_p_pearson <- NA
  w_mad <- NA
  w_rmsd <- NA
  w_spearman <- NA
  w_y_pearson <- NA
  gridres <- NULL
  if (grid.search) {
    message("Grid search for ENSEMBLE weight ...")
    gridmat <- getSearchGrid(lengthby = search.length, nparam = length(prop.list))
    ptm <- proc.time()
    search.prop <- NULL
    if (!is.null(truep)) {
      message("Searching according to proportion--Pearson Correlation...")
      search.prop <- t(apply(gridmat, 1, function(x) {
        temp <- wt_prop(x, proplist = combo)
        w_eval <- SCDC_peval(ptrue = as.matrix(truep), 
                             pest = temp, pest.names = c("SCDC"), select.ct = ct.sub)
        w_eval$evals.table
      }))
      colnames(search.prop) <- c("RMSD", "mAD", "R")
      w_p_pearson <- gridmat[which.max(search.prop[, 3]), 
                             ]
    }
    message("Searching according to bulk expression measurement...")
    search.y <- t(apply(gridmat, 1, function(x) {
      temp <- wt_y(x, y.list = y.list)
      c(sqrt(mean((yv - temp)^2)), mean(abs(yv - temp)), 
        cor(yv, temp, method = "pearson"), cor(yv, temp, 
                                               method = "spearman"))
    }))
    colnames(search.y) <- c("RMSD_Y", "mAD_Y", "Pearson_Y", 
                            "Spearman_Y")
    ptm2 <- proc.time() - ptm
    message("Grid search used", ptm2[3], " seconds.")
    w_rmsd <- gridmat[which.min(search.y[, 1]), ]
    w_mad <- gridmat[which.min(search.y[, 2]), ]
    w_y_pearson <- gridmat[which.max(search.y[, 3]), ]
    w_spearman <- gridmat[which.max(search.y[, 4]), ]
  }
  if (!is.null(truep)) {
    weight.mat <- rbind(sse.wt, sae.wt, rmsd.wt, w_lad, w_nnls, 
                        w_p_pearson, w_mad, w_rmsd, w_spearman)
    rownames(weight.mat) <- c("inverse SSE", "inverse SAE", 
                              "inverse RMSD", "LAD", "NNLS", "Pearson_prop", "mAD_Y", 
                              "RMSD_Y", "Spearman_Y")
    eval.prop <- t(apply(weight.mat, 1, function(x) {
      temp <- wt_prop(x, proplist = combo)
      w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp, 
                           pest.names = c("SCDC"), select.ct = ct.sub)
      w_eval$evals.table
    }))
    colnames(eval.prop) <- c("RMSD_prop", "mAD_prop", "R_prop")
    eval.y <- t(apply(weight.mat, 1, function(x) {
      temp <- wt_y(x, y.list = y.list)
      c(sqrt(mean((yv - temp)^2)), mean(abs(yv - temp)), 
        cor(yv, temp, method = "pearson"), cor(yv, temp, 
                                               method = "spearman"))
    }))
    colnames(eval.y) <- c("RMSD_Y", "mAD_Y", "Pearson_Y", 
                          "Spearman_Y")
    if (grid.search) {
      gridres <- cbind(search.prop, search.y, gridmat)
    }
    out <- cbind(round(weight.mat, 2), eval.y, eval.prop)
  }
  else {
    weight.mat <- rbind(sse.wt, sae.wt, rmsd.wt, w_lad, w_nnls, 
                        w_mad, w_rmsd, w_spearman)
    rownames(weight.mat) <- c("inverse SSE", "inverse SAE", 
                              "inverse RMSD", "LAD", "NNLS", "mAD_Y", "RMSD_Y", 
                              "Spearman_Y")
    eval.y <- t(apply(weight.mat, 1, function(x) {
      temp <- wt_y(x, y.list = y.list)
      c(sqrt(mean((yv - temp)^2)), mean(abs(yv - temp)), 
        cor(yv, temp, method = "pearson"), cor(yv, temp, 
                                               method = "spearman"))
    }))
    colnames(eval.y) <- c("RMSD_Y", "mAD_Y", "Pearson_Y", 
                          "Spearman_Y")
    if (grid.search) {
      gridres <- cbind(search.y, gridmat)
    }
    out <- cbind(round(weight.mat, 2), eval.y)
  }
  return(list(w_table = out, prop.list = prop.list, prop.only = combo, 
              gridres = gridres))
}

getSearchGrid <- function (lengthby, nparam) {
  wlist <- list()
  for (i in 1:nparam) {
    wlist[[i]] <- seq(0, 1 + lengthby, by = lengthby)
  }
  w1grid <- round(expand.grid(wlist), digits = 2)
  w1grid <- w1grid[round(rowSums(w1grid), digits = 2) == 1, 
                   ]
  colnames(w1grid) <- paste("w", 1:nparam, sep = "")
  return(w1grid)
}

wt_extP <- function(wt,proplist){
  wt <- as.numeric(wt)
  combo.list <- list()
  for (one in proplist) {
    combo.list[[length(combo.list)+1]] <- one$pV# * wt[length(combo.list)+1]
  }
  combo.pV <- apply(simplify2array(combo.list),1:2,min)
  return(combo.pV)
  
}
wt_prop <- function (wt, proplist){
  wt <- as.numeric(wt)
  combo.list <- list()
  for (i in 1:length(proplist)) {
    combo.list[[i]] <- proplist[[i]] * wt[i]
  }
  combo.prop <- Reduce("+", combo.list)
  return(combo.prop)
}
wt_y <- function (wt, y.list = y.list) {
  wt <- as.numeric(wt)
  combo.list <- list()
  for (i in 1:ncol(y.list)) {
    combo.list[[i]] <- y.list[, i] * wt[i]
  }
  combo.y <- Reduce("+", combo.list)
  return(combo.y)
}
wt_cover <- function(proplist,scG){
  geneL <- c()
  for(i in proplist) geneL <- unique(c(geneL,rownames(i$yhat)))
  return(round(100*sum(geneL%in%scG)/length(scG),1))
  
}
wt_pval <- function(A,b,beta){
  Init <- setNames(beta,paste("b",1:length(beta),sep=""))
  formu <- paste("bulk~",paste(paste(names(Init),"*",colnames(A),sep=""),collapse="+"),sep="")
  Data <- as.data.frame(cbind(A,bulk=as.vector(b)))
  fit <- try(nls(as.formula(formu),Data,Init,algorithm="port",lower=rep(0,length(Init))),silent=T)
  pV <- setNames(rep(0.0499,length(beta)+1),c(colnames(A),"overallP"))
  if("convergence"%in%names(fit)){
    a <- summary(fit)
    pV <- setNames(c(a$coefficients[,"Pr(>|t|)"],cor.test(predict(fit),b)$p.value),c(colnames(A),"overallP"))
  }else{
    print(fit)
  }
  return(pV)
}

SCDC_peval <- function (ptrue, pest, pest.names, select.ct = NULL) {
  if (!is.list(pest)) {
    pest <- list(pest)
  }
  if (!is.data.frame(ptrue)) {
    ptrue <- as.data.frame.matrix(ptrue)
  }
  n_est <- length(pest)
  sample_names <- lapply(pest, rownames)
  ctype_names <- lapply(pest, colnames)
  sample_common <- Reduce(intersect, sample_names)
  ctype_common <- Reduce(intersect, ctype_names)
  celltype <- intersect(colnames(ptrue), ctype_common)
  if (!is.null(select.ct)) {
    celltype <- intersect(celltype, select.ct)
  }
  sample <- intersect(rownames(ptrue), sample_common)
  N <- length(sample)
  K <- length(celltype)
  if (N < 1) {
    stop("No common Subjects! Check rowname!")
  }
  if (K <= 1) {
    stop("Not enough cell types!")
  }
  ptrue.use <- ptrue[intersect(rownames(ptrue), sample), intersect(colnames(ptrue), 
                                                                   celltype)]
  ptrue.use <- as.data.frame.matrix(ptrue.use/apply(ptrue.use, 
                                                    1, sum))
  ptrue.use[is.na(ptrue.use)] <- 0
  evals <- lapply(pest, function(xx) {
    pest.use <- xx[intersect(rownames(xx), sample), intersect(colnames(xx), 
                                                              celltype)]
    pest.use <- as.data.frame.matrix(pest.use/apply(pest.use, 
                                                    1, sum))
    pest.use <- pest.use[rownames(ptrue.use), colnames(ptrue.use)]
    RMSD_bysample <- round(sqrt(rowMeans((ptrue.use - pest.use)^2)), 
                           digits = 5)
    mAD_bysample <- round(rowMeans(abs(ptrue.use - pest.use)), 
                          digits = 5)
    Pearson_bysample <- sapply(1:nrow(ptrue.use), function(ss) {
      round(cor(c(as.matrix(ptrue.use[ss, ])), c(as.matrix(pest.use[ss, 
                                                                    ]))), digits = 5)
    })
    RMSD <- round(sqrt(mean(as.matrix((ptrue.use - pest.use)^2), 
                            na.rm = T)), digits = 5)
    mAD <- round(mean(as.matrix(abs(ptrue.use - pest.use)), 
                      na.rm = T), digits = 5)
    Pearson <- round(cor(c(as.matrix(ptrue.use)), c(as.matrix(pest.use))), 
                     digits = 4)
    return(list(pest.use = pest.use, RMSD_bysample = RMSD_bysample, 
                mAD_bysample = mAD_bysample, Pearson_bysample = Pearson_bysample, 
                RMSD = RMSD, mAD = mAD, Pearson = Pearson))
  })
  evals.table <- NULL
  for (l in 1:n_est) {
    evals.table <- rbind(evals.table, c(evals[[l]]$RMSD, 
                                        evals[[l]]$mAD, evals[[l]]$Pearson))
  }
  colnames(evals.table) <- c("RMSD", "mAD", "R")
  rownames(evals.table) <- pest.names
  pearson.sample.table <- NULL
  for (l in 1:n_est) {
    pearson.sample.table <- rbind(pearson.sample.table, evals[[l]]$Pearson_bysample)
  }
  rownames(pearson.sample.table) <- pest.names
  colnames(pearson.sample.table) <- rownames(ptrue.use)
  RMSD.sample.table <- NULL
  for (l in 1:n_est) {
    RMSD.sample.table <- rbind(RMSD.sample.table, evals[[l]]$RMSD_bysample)
  }
  rownames(RMSD.sample.table) <- pest.names
  colnames(RMSD.sample.table) <- rownames(ptrue.use)
  mAD.sample.table <- NULL
  for (l in 1:n_est) {
    mAD.sample.table <- rbind(mAD.sample.table, evals[[l]]$mAD_bysample)
  }
  rownames(mAD.sample.table) <- pest.names
  colnames(mAD.sample.table) <- rownames(ptrue.use)
  return(list(evals = evals, evals.table = evals.table, pearson.sample.table = pearson.sample.table, 
              RMSD.sample.table = RMSD.sample.table, mAD.sample.table = mAD.sample.table))
}

SCDC_parrel_deconv <- function(bulkID){
	
	return(bplapply(bulkID,function(i){#bplapply#sapply
		nu <- 1e-4
		epsilon <- 0.01
		iter.max <- 1000
		ALS.S <- SCDCdeconvObj$ALS.S

		basis.mvw.temp <- SCDCdeconvObj$basis.mvw
		xbulk.temp <- SCDCdeconvObj$xbulk[, i] * 100
		sigma.temp <- SCDCdeconvObj$sigma
		#message(paste(i, "has common genes",sum(SCDCdeconvObj$xbulk[, i] != 0), "..."))
		lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
		delta <- lm$residuals
		wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 * 
											 	t(sigma.temp)))
		x.wt <- xbulk.temp * sqrt(wt.gene)
		b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
		lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
		prop.wt <- lm.wt$x/sum(lm.wt$x)
		delta <- lm.wt$residuals
		for (iter in 1:iter.max) {
			wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * 
												  	ALS.S)^2 * t(sigma.temp)))
			x.wt <- xbulk.temp * sqrt(wt.gene)
			b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), 
						  "*")
			lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
			delta.new <- lm.wt$residuals
			prop.wt.new <- lm.wt$x/sum(lm.wt$x)
			if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
				prop.wt <- prop.wt.new
				delta <- delta.new
				#R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% 
				#                as.matrix(lm.wt$x))/var(xbulk.temp)
				#message("WNNLS Converged at iteration ", iter)
				break
			}
			prop.wt <- prop.wt.new
			delta <- delta.new
		}
		pV.tmp <- wt_pval(b.wt,x.wt,lm.wt$x)## adding to provide pvalues!
		R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
		yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
		
		return(list(R2=R2,prop.wt=prop.wt,yhat.temp=yhat.temp,pV.tmp=pV.tmp))
		
	}))
}

SCDC_prop <- function (bulk.eset, sc.eset, ct.varname, sample, ct.sub, iter.max = 1000, 
                       nu = 1e-04, epsilon = 0.01, truep = NULL, weight.basis = T, 
                       Transform_bisque = F, ...) {
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[,ct.varname]))
  sc.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, 
                         sample = sample)
  commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  if (weight.basis) {
    basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
  }
  else {
    basis.mvw <- sc.basis$basis[commongenes, ct.sub]
  }
  scdcP <- NULL## adding to provide pvalues!
  if (Transform_bisque) {
    GenerateSCReference <- function(sc.eset, ct.sub) {
      cell.labels <- base::factor(sc.eset[[ct.sub]])
      all.cell.types <- base::levels(cell.labels)
      aggr.fn <- function(ct.sub) {
        base::rowMeans(Biobase::exprs(sc.eset)[, cell.labels == 
                                                 ct.sub, drop = F])
      }
      template <- base::numeric(base::nrow(sc.eset))
      sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
      return(sc.ref)
    }
    sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes, 
                                                       , drop = F]
    ncount <- table(sc.eset@phenoData@data[, sample], sc.eset@phenoData@data[, 
                                                                             ct.varname])
    true.prop <- ncount/rowSums(ncount, na.rm = T)
    sc.props <- round(true.prop[complete.cases(true.prop), 
                                ], 2)
    Y.train <- sc.ref %*% t(sc.props[, colnames(sc.ref)])
    dim(Y.train)
    X.pred <- exprs(bulk.eset)[commongenes, ]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    SemisupervisedTransformBulk <- function(gene, Y.train, 
                                            X.pred) {
      Y.train.scaled <- base::scale(Y.train[gene, , drop = T])
      Y.center <- base::attr(Y.train.scaled, "scaled:center")
      Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
      n <- base::length(Y.train.scaled)
      shrink.scale <- base::sqrt(base::sum((Y.train[gene, 
                                                    , drop = T] - Y.center)^2)/n + 1)
      X.pred.scaled <- base::scale(X.pred[gene, , drop = T])
      Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + 
                               Y.center, dimnames = base::list(base::colnames(X.pred), 
                                                               gene))
      return(Y.pred)
    }
    Y.pred <- base::matrix(base::vapply(X = commongenes, 
                                        FUN = SemisupervisedTransformBulk, FUN.VALUE = template, 
                                        Y.train, X.pred, USE.NAMES = TRUE), nrow = base::length(sample.names))
    indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
      base::anyNA(column)
    })
    if (base::any(indices)) {
      if (sum(!indices) == 0) {
        base::stop("Zero genes left for decomposition.")
      }
      Y.pred <- Y.pred[, !indices, drop = F]
      sc.ref <- sc.ref[!indices, , drop = F]
    }
    results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
      sol <- lsei::pnnls(sc.ref, b, sum = 1)
      return(sol$x)
    }))
    prop.est.mvw <- t(results)
    colnames(prop.est.mvw) <- colnames(sc.ref)
    rownames(prop.est.mvw) <- colnames(bulk.eset)
    yhat <- sc.ref %*% results
    colnames(yhat) <- colnames(bulk.eset)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }
  else {
    xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
    sigma <- sc.basis$sigma[commongenes, ct.sub]
    ALS.S <- sc.basis$sum.mat[ct.sub]
    N.bulk <- ncol(bulk.eset)
    valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) ==0) & (!is.na(ALS.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    basis.mvw <- basis.mvw[, valid.ct]
    ALS.S <- ALS.S[valid.ct]
    sigma <- sigma[, valid.ct]
    prop.est.mvw <- NULL
    yhat <- NULL
    yhatgene.temp <- rownames(basis.mvw)
    if(T){
      if(F){
      allRes <- bplapply(colnames(xbulk),function(i){#bplapply#sapply
        basis.mvw.temp <- basis.mvw
        xbulk.temp <- xbulk[, i] * 100
        sigma.temp <- sigma
        #message(paste(i, "has common genes",sum(xbulk[, i] != 0), "..."))
        lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
        delta <- lm$residuals
        wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 * 
                                               t(sigma.temp)))
        x.wt <- xbulk.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        prop.wt <- lm.wt$x/sum(lm.wt$x)
        delta <- lm.wt$residuals
        for (iter in 1:iter.max) {
          wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * 
                                                  ALS.S)^2 * t(sigma.temp)))
          x.wt <- xbulk.temp * sqrt(wt.gene)
          b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), 
                        "*")
          lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
          delta.new <- lm.wt$residuals
          prop.wt.new <- lm.wt$x/sum(lm.wt$x)
          if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
            prop.wt <- prop.wt.new
            delta <- delta.new
            #R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% 
            #                as.matrix(lm.wt$x))/var(xbulk.temp)
            #message("WNNLS Converged at iteration ", iter)
            break
          }
          prop.wt <- prop.wt.new
          delta <- delta.new
        }
        pV.tmp <- wt_pval(b.wt,x.wt,lm.wt$x)## adding to provide pvalues!
        R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
        yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
        
        return(list(R2=R2,prop.wt=prop.wt,yhat.temp=yhat.temp,pV.tmp=pV.tmp))
        
      })#,BPPARAM=MulticoreParam(workers=min(multicoreWorkers(),32),tasks=1)
      }
      SCDCdeconvObj <<- list(xbulk=xbulk,basis.mvw=basis.mvw,sigma=sigma,ALS.S=ALS.S)
      #message("SCDC_prop start1")
      allRes <- SCDC_parrel_deconv(colnames(xbulk))
      rm(SCDCdeconvObj,pos=".GlobalEnv")
      #message("SCDC_prop end1")
      for(one in allRes){
        prop.est.mvw <- rbind(prop.est.mvw, one$prop.wt)
        yhatgene.temp <- intersect(rownames(one$yhat.temp), yhatgene.temp)
        yhat <- cbind(yhat[yhatgene.temp, ],one$yhat.temp[yhatgene.temp,])
        scdcP <- cbind(scdcP,one$pV.tmp)## adding to provide pvalues!
      }
    }else{
      for (i in 1:N.bulk) {
        basis.mvw.temp <- basis.mvw
        xbulk.temp <- xbulk[, i] * 100
        sigma.temp <- sigma
        message(paste(colnames(xbulk)[i], "has common genes", 
                      sum(xbulk[, i] != 0), "..."))
        lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
        delta <- lm$residuals
        wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 * 
                                               t(sigma.temp)))
        x.wt <- xbulk.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        prop.wt <- lm.wt$x/sum(lm.wt$x)
        delta <- lm.wt$residuals
        for (iter in 1:iter.max) {
          wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * 
                                                  ALS.S)^2 * t(sigma.temp)))
          x.wt <- xbulk.temp * sqrt(wt.gene)
          b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), 
                        "*")
          lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
          delta.new <- lm.wt$residuals
          prop.wt.new <- lm.wt$x/sum(lm.wt$x)
          if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
            prop.wt <- prop.wt.new
            delta <- delta.new
            R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% 
                            as.matrix(lm.wt$x))/var(xbulk.temp)
            message("WNNLS Converged at iteration ", iter)
            break
          }
          prop.wt <- prop.wt.new
          delta <- delta.new
        }
        pV.tmp <- wt_pval(b.wt,x.wt,lm.wt$x)## adding to provide pvalues!
        
        R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
        prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
        yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
        yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
        yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp,])
        scdcP <- cbind(scdcP,pV.tmp)## adding to provide pvalues!
        #browser()
      }
    }
    #browser()
    dimnames(prop.est.mvw) <- list(colnames(xbulk),colnames(basis.mvw))
    colnames(yhat) <- colnames(scdcP) <- colnames(xbulk)## adding to provide pvalues!
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw, 
              yhat = yhat, yeval = yeval, peval = peval, pV=scdcP))## adding to provide pvalues!
}

SCDC_prop_ONE <- function (bulk.eset, sc.eset, ct.varname, sample, truep = NULL, 
                           ct.sub, iter.max = 2000, nu = 1e-10, epsilon = 0.01, weight.basis = T, 
                           ...) {
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, 
                             ct.varname = ct.varname, sample = sample)
  if (weight.basis) {
    basis <- sc.basis$basis.mvw
  }
  else {
    basis <- sc.basis$basis
  }
  commongenes <- intersect(rownames(basis), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  basis.mvw <- basis[commongenes, ct.sub]
  xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  prop.est.mvw <- NULL
  yhat <- NULL
  yhatgene.temp <- rownames(basis.mvw)
  scdcP <- NULL## adding to provide pvalues!
  for (i in 1:N.bulk) {
    xbulk.temp <- xbulk[, i]
    message(paste(colnames(xbulk)[i], "has common genes", 
                  sum(xbulk[, i] != 0), "..."))
    lm <- nnls::nnls(A = basis.mvw, b = xbulk.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp * sqrt(wt.gene)
    b.wt <- sweep(basis.mvw, 1, sqrt(wt.gene), "*")
    lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals
    for (iter in 1:iter.max) {
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)
      if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("WNNLS Converged at iteration ", iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }
    scdcP <- cbind(scdcP,wt_pval(b.wt,x.wt,lm.wt$x))## adding to provide pvalues!
    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
    yhat.temp <- basis.mvw %*% as.matrix(lm.wt$x)
    yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
    yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp, 
                                                   ])
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(bulk.eset)
  colnames(yhat) <- colnames(scdcP) <- colnames(bulk.eset) ## adding to provide pvalues!
  yobs <- exprs(bulk.eset)
  yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
  peval <- NULL
  if (!is.null(truep)) {
    if (all(rownames(truep) == rownames(prop.est.mvw))) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
    else {
      message("Your input sample names for proportion matrix and bulk.eset do not match! Please make sure sample names match.")
    }
  }
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw, 
              yhat = yhat, yeval = yeval, peval = peval,pV=scdcP))## adding to provide pvalues!
}
SCDC_basis_ONE <- function (x, ct.sub = NULL, ct.varname, sample) {
  if (is.null(ct.sub)) {
    ct.sub <- unique(x@phenoData@data[, ct.varname])[!is.na(unique(x@phenoData@data[, 
                                                                                    ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[, x@phenoData@data[, ct.varname] %in% ct.sub]
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0, ]
  countmat <- exprs(x.sub)
  ct.id <- x.sub@phenoData@data[, ct.varname]
  sample.id <- x.sub@phenoData@data[, sample]
  ct_sample.id <- paste(ct.id, sample.id, sep = "%")
  mean.mat <- sapply(unique(ct_sample.id), function(id) {
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y, 1, sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call("rbind", strsplit(unique(ct_sample.id), 
                                       split = "%"))
  sum.mat2 <- sapply(unique(sample.id), function(sid) {
    sapply(unique(ct.id), function(id) {
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% 
                               sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <- rowMeans(sum.mat2, na.rm = T)
  basis <- sapply(unique(mean.id[, 1]), function(id) {
    z <- sum.mat[mean.id[, 1]]
    mean.mat.z <- t(t(mean.mat) * z)
    y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
    apply(y, 1, mean, na.rm = TRUE)
  })
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, na.rm = TRUE)
    y/median(y, na.rm = T)
  }
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid, 
                   drop = FALSE]
      apply(y, 1, var, na.rm = T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)
  q15 <- apply(var.adj, 2, quantile, probs = 0.15, na.rm = T)
  q85 <- apply(var.adj, 2, quantile, probs = 0.85, na.rm = T)
  var.adj.q <- as.matrix(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)
  }) + 1e-04)
  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id) {
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q), "/")
    apply(yy, 1, sum, na.rm = TRUE)/sum(yy)
  })
  basis.mvw <- sapply(unique(mean.id[, 1]), function(id) {
    z <- sum.mat[mean.id[, 1]]
    mean.mat.z <- t(t(mean.mat.mvw) * z)
    y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
    apply(y, 1, mean, na.rm = TRUE)
  })
  basis.mvw <- basis.mvw[, ct.sub]
  sigma <- NULL
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]
  return(list(basis = basis, sum.mat = sum.mat, sigma = sigma, 
              basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

SCDC_basis <- function (x, ct.sub = NULL, ct.varname, sample) {
  bePar <- T
  if (is.null(ct.sub)) {
    ct.sub <- unique(x@phenoData@data[, ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[, x@phenoData@data[, ct.varname] %in% ct.sub]
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0, ]
  if(bePar){
  	countmat <<- exprs(x.sub)
  	ct.id <<- droplevels(as.factor(x.sub@phenoData@data[, ct.varname]))
  	sample.id <<- as.character(x.sub@phenoData@data[, sample])
  	ct_sample.id <<- paste(ct.id, sample.id, sep = "%")
  	rm(x,x.sub)
  	#message("start1")
    tmp <- bplapply(unique(ct_sample.id), function(id) {
      #message(id)
      y = as.matrix(countmat[, ct_sample.id %in% id,drop=F])
      return(apply(y, 1, sum, na.rm = TRUE)/sum(y))
    })
    #message("end1")
    mean.mat <<- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(mean.mat) <<- unique(ct_sample.id)
  }else{
  	countmat <- exprs(x.sub)
  	ct.id <- droplevels(as.factor(x.sub@phenoData@data[, ct.varname]))
  	sample.id <- as.character(x.sub@phenoData@data[, sample])
  	ct_sample.id <- paste(ct.id, sample.id, sep = "%")
  	mean.mat <- sapply(unique(ct_sample.id), function(id) {
      y = as.matrix(countmat[, ct_sample.id %in% id,drop=F])
      apply(y, 1, sum, na.rm = TRUE)/sum(y)
    })
  }
  mean.id <- do.call("rbind", strsplit(unique(ct_sample.id),split = "%"))
  if(bePar){
  	#message("start2")
    tmp <- bplapply(unique(mean.id[, 1]), function(id) {
      y = mean.mat[, mean.id[, 1] %in% id,drop=F]
      return(apply(y, 1, var, na.rm = TRUE))
    })
    #message("end2")
    sigma <- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(sigma) <- unique(mean.id[, 1])
  }else{
    sigma <- sapply(unique(mean.id[, 1]), function(id) {
      y = mean.mat[, mean.id[, 1] %in% id,drop=F]
      apply(y, 1, var, na.rm = TRUE)
    })
  }
  
  if(bePar){
  	#message("start3")
    tmp <- bplapply(unique(sample.id), function(sid) {
      sapply(unique(ct.id), function(id) {
        y = as.matrix(countmat[, ct.id %in% id & sample.id %in% 
                                 sid,drop=F])
        return(sum(y)/ncol(y))
      })
    })
    #message("end3")
    sum.mat2 <- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(sum.mat2) <- unique(sample.id)
  }else{
    sum.mat2 <- sapply(unique(sample.id), function(sid) {
      sapply(unique(ct.id), function(id) {
        y = as.matrix(countmat[, ct.id %in% id & sample.id %in% 
                                 sid,drop=F])
        sum(y)/ncol(y)
      })
    })
  }  
  
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <<- rowMeans(sum.mat2, na.rm = T)
  
  if(bePar){
  	#message("start4")
    tmp <- bplapply(unique(mean.id[, 1]), function(id) {
      z <- sum.mat[mean.id[, 1]]
      mean.mat.z <- t(t(mean.mat) * z)
      y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id,drop=F])
      return(apply(y, 1, mean, na.rm = TRUE))
    })
    #message("end4")
    basis <- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(basis) <- unique(mean.id[, 1])
  }else{
    basis <- sapply(unique(mean.id[, 1]), function(id) {
      z <- sum.mat[mean.id[, 1]]
      mean.mat.z <- t(t(mean.mat) * z)
      y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id,drop=F])
      apply(y, 1, mean, na.rm = TRUE)
    })
  }  
  
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, ...)
    yM <- median(y,...)
    if(yM==0)yM <- min(y[y>0])
    y/yM
  }
  if(bePar){
  	#message("start5")
    tmp <- bplapply(unique(sample.id), function(sid) {
      return(my.max(sapply(unique(ct.id), function(id) {
        y = countmat[, ct.id %in% id & sample.id %in% sid, 
                     drop = FALSE]
        return(apply(y, 1, var, na.rm = T))
      }), na.rm = T))
    })
    #message("end5")
    var.adj <- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(var.adj) <- unique(sample.id)
  }else{
    var.adj <- sapply(unique(sample.id), function(sid) {
      my.max(sapply(unique(ct.id), function(id) {
        y = countmat[, ct.id %in% id & sample.id %in% sid, 
                     drop = FALSE]
        apply(y, 1, var, na.rm = T)
      }), na.rm = T)
    })
  }  
  
  colnames(var.adj) <- unique(sample.id)
  q15 <- apply(var.adj, 2, quantile, probs = 0.15, na.rm = T)
  q85 <- apply(var.adj, 2, quantile, probs = 0.85, na.rm = T)
  
  var.adj.q <- t(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)
  })) + 1e-04
  message("Creating Basis Matrix adjusted for maximal variance weight")
  if(bePar){
  	#message("start6")
    tmp <- bplapply(unique(ct_sample.id), function(id) {
      sid = unlist(strsplit(id, "%"))[2]
      y = as.matrix(countmat[, ct_sample.id %in% id])
      yy = sweep(y, 1, sqrt(var.adj.q[, sid]), "/")
      return(apply(yy, 1, sum, na.rm = TRUE)/sum(yy))
    })
    #message("end6")
    mean.mat.mvw <<- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(mean.mat.mvw) <<- unique(ct_sample.id)
  }else{
    mean.mat.mvw <- sapply(unique(ct_sample.id), function(id) {
      sid = unlist(strsplit(id, "%"))[2]
      y = as.matrix(countmat[, ct_sample.id %in% id])
      yy = sweep(y, 1, sqrt(var.adj.q[, sid]), "/")
      apply(yy, 1, sum, na.rm = TRUE)/sum(yy)
    })
  }   
  if(bePar){
  	#message("start7")
    tmp <- bplapply(unique(mean.id[, 1]), function(id) {
      z <- sum.mat[mean.id[, 1]]
      mean.mat.z <- t(t(mean.mat.mvw) * z)
      y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
      return(apply(y, 1, mean, na.rm = TRUE))
    })
    #message("end7")
    basis.mvw <- sapply(tmp,function(one)return(one))
    rm(tmp)
    colnames(basis.mvw) <- unique(mean.id[, 1])
	rm(countmat,ct.id,sample.id,ct_sample.id,mean.mat,mean.mat.mvw,pos=".GlobalEnv")
	gc()    
  }else{
    basis.mvw <- sapply(unique(mean.id[, 1]), function(id) {
      z <- sum.mat[mean.id[, 1]]
      mean.mat.z <- t(t(mean.mat.mvw) * z)
      y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
      apply(y, 1, mean, na.rm = TRUE)
    })
  }   
  
  basis.mvw <- basis.mvw[, ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat, sigma = sigma, 
              basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

getCPM0 <- function (x, verbose = F) {
  if (is.null(dim(x))) {
    if (verbose) {
      message("Normalizing a vector instead of a matrix")
    }
    vec = as.matrix(x/sum(x))
    vec
  }
  else {
    cpm <- t(t(x)/apply(x, 2, sum))
    cpm
  }
}

SCDC_yeval <- function (y, yest, yest.names = NULL) {
  if (!is.list(yest)) {
    yest <- list(yest)
  }
  n_est <- length(yest)
  evals <- lapply(yest, function(xx) {
    if (!is.null(xx)) {
      if (dim(xx)[2] > 1) {
        g.use = intersect(rownames(y), rownames(xx))
        y.norm <- getCPM0(y[g.use, ])
        x = xx[g.use, colnames(y)]
        yest.norm = getCPM0(x)
        spearmany <- round(cor(c(yest.norm), c(y.norm), 
                               method = "spearman"), digits = 5)
        RMSDy = round(sqrt(mean((yest.norm - y.norm)^2)), 
                      digits = 7)
        mADy = round(mean(abs(yest.norm - y.norm)), digits = 8)
        spearmany_bysample <- sapply(1:ncol(x), function(ss) {
          round(cor(yest.norm[, ss], y.norm[, ss], method = "spearman"), 
                digits = 4)
        })
        RMSDy_bysample = round(sqrt(colMeans((yest.norm[, 
                                                        ] - y.norm[, ])^2)), digits = 6)
        mADy_bysample = round(colMeans(abs(yest.norm - 
                                             y.norm)), digits = 6)
      }
      else {
        g.use = intersect(rownames(y), rownames(xx))
        y.norm <- getCPM0(y[g.use, ])
        x = xx[g.use, colnames(y)]
        yest.norm = getCPM0(x)
        spearmany <- round(cor(c(yest.norm), c(y.norm), 
                               method = "spearman"), digits = 5)
        RMSDy = round(sqrt(mean((yest.norm - y.norm)^2)), 
                      digits = 7)
        mADy = round(mean(abs(yest.norm - y.norm)), digits = 8)
        spearmany_bysample <- spearmany
        RMSDy_bysample = RMSDy
        mADy_bysample = mADy
      }
    }
    else if (is.null(xx)) {
      spearmany = NA
      RMSDy = NA
      mADy = NA
      spearmany_bysample = NA
      RMSDy_bysample = NA
      mADy_bysample = NA
    }
    return(list(spearmany = spearmany, RMSDy = RMSDy, mADy = mADy, 
                spearmany_bysample = spearmany_bysample, RMSDy_bysample = RMSDy_bysample, 
                mADy_bysample = mADy_bysample))
  })
  yevals.table <- NULL
  for (l in 1:n_est) {
    yevals.table <- rbind(yevals.table, c(evals[[l]]$spearmany, 
                                          evals[[l]]$RMSDy, evals[[l]]$mADy))
  }
  colnames(yevals.table) <- c("spearman_Y", "RMSD_Y", "mAD_Y")
  rownames(yevals.table) <- yest.names
  spearmany.sample.table <- NULL
  for (l in 1:n_est) {
    spearmany.sample.table <- rbind(spearmany.sample.table, 
                                    evals[[l]]$spearmany_bysample)
  }
  rownames(spearmany.sample.table) <- yest.names
  colnames(spearmany.sample.table) <- colnames(y)
  RMSDy.sample.table <- NULL
  for (l in 1:n_est) {
    RMSDy.sample.table <- rbind(RMSDy.sample.table, evals[[l]]$RMSDy_bysample)
  }
  rownames(RMSDy.sample.table) <- yest.names
  colnames(RMSDy.sample.table) <- colnames(y)
  mADy.sample.table <- NULL
  for (l in 1:n_est) {
    mADy.sample.table <- rbind(mADy.sample.table, evals[[l]]$mADy_bysample)
  }
  rownames(mADy.sample.table) <- yest.names
  colnames(mADy.sample.table) <- colnames(y)
  return(list(evals = evals, yevals.table = yevals.table, spearmany.sample.table = spearmany.sample.table, 
              RMSDy.sample.table = RMSDy.sample.table, mADy.sample.table = mADy.sample.table))
}

