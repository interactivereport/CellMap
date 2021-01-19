##########################
## deconv.Bisque.R
##
##########################
loadBisque <- function(){
  if(!require(Biobase)){
    BiocManager::install("Biobase")
    if(!require(Biobase)) stop("Cannot install Biobase!")
  }
  if(!require(lsei)){
    install.packages("lsei")
    if(!require(lsei)) stop("Cannot install lsei!")
  }
}
suppressMessages(loadBisque())

deconv.Bisque <- function(bulk,feature){
  SC <- feature$SC
  rm(feature)
  comp <- ReferenceBasedDecomposition(ExpressionSet(as.matrix(bulk)), SC,cell.types='cellType',subject.names='dID', use.overlap=F)
  pV <- comp$bulk.props
  pV[,] <- 0.049
  overallP <- setNames(rep(0.049,ncol(pV)),colnames(pV))
  return(list(composition=comp$bulk.props,
              compoP=pV,
              overallP=overallP,
              coverR=round(100*length(comp$genes.used)/nrow(SC),1),
              rawComp=NULL,
              rawSets=NULL,
              missingF=""))
}

ReferenceBasedDecomposition <- function (bulk.eset, sc.eset, markers = NULL, cell.types = "cellType", 
          subject.names = "SubjectName", use.overlap = TRUE, verbose = TRUE) 
{
  if ((!methods::is(sc.eset, "ExpressionSet")) || (!methods::is(bulk.eset, 
                                                                "ExpressionSet"))) {
    base::stop("Expression data should be in ExpressionSet")
  }
  else if (!cell.types %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Cell type label \"%s\" ", cell.types), 
               "not found in single-cell ExpressionSet varLabels.")
  }
  else if (!subject.names %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Individual label \"%s\"", subject.names), 
               " not found in single-cell ExpressionSet varLabels.")
  }
  n.sc.individuals <- base::length(base::levels(base::factor(sc.eset[[subject.names]])))
  if (n.sc.individuals == 1) {
    base::stop("Only one individual detected in single-cell data. At least ", 
               "two subjects are needed (three or more recommended).")
  }
  else if (n.sc.individuals == 2) {
    base::warning("Only two individuals detected in single-cell data. While ", 
                  "Bisque will run, we recommend at least three subjects for", 
                  " reliable performance.")
  }
  n.cell.types <- base::length(base::levels(base::factor(sc.eset[[cell.types]])))
  if (n.cell.types == 1) {
    base::stop("Single-cell pheno data indicates only one cell type", 
               " present. No need for decomposition.")
  }
  if (verbose) {
    base::message(base::sprintf("Decomposing into %i cell types.", 
                                n.cell.types))
  }
  if (use.overlap) {
    samples <- GetOverlappingSamples(sc.eset, bulk.eset, 
                                     subject.names, verbose)
  }
  if (base::is.null(markers)) {
    markers <- Biobase::featureNames(sc.eset)
  }
  else {
    markers <- base::unique(base::unlist(markers))
  }
  genes <- GetOverlappingGenes(sc.eset, bulk.eset, markers, 
                               verbose)
  sc.eset <- Biobase::ExpressionSet(assayData = Biobase::exprs(sc.eset)[genes, 
                                                                        ], phenoData = sc.eset@phenoData)
  bulk.eset <- Biobase::ExpressionSet(assayData = Biobase::exprs(bulk.eset)[genes, 
                                                                            ], phenoData = bulk.eset@phenoData)
  if (verbose) {
    base::message("Converting single-cell counts to CPM and ", 
                  "filtering zero variance genes.")
  }
  sc.eset <- CountsToCPM(sc.eset)
  sc.eset <- FilterZeroVarianceGenes(sc.eset, verbose)
  if (verbose) {
    base::message("Converting bulk counts to CPM and filtering", 
                  " unexpressed genes.")
  }
  bulk.eset <- CountsToCPM(bulk.eset)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  genes <- base::intersect(Biobase::featureNames(sc.eset), 
                           Biobase::featureNames(bulk.eset))
  if (base::length(genes) == 0) {
    base::stop("Zero genes remaining after filtering and ", 
               "intersecting bulk, single-cell, and marker genes.")
  }
  if (verbose) {
    n.cells <- base::ncol(sc.eset)
    base::message("Generating single-cell based reference from ", 
                  sprintf("%i cells.\n", n.cells))
  }
  sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes, 
                                                     , drop = F]
  sc.props <- CalculateSCCellProportions(sc.eset, subject.names, 
                                         cell.types)
  sc.props <- sc.props[base::colnames(sc.ref), , drop = F]
  if (use.overlap) {
    if (verbose) {
      base::message("Learning bulk transformation from overlapping samples.")
    }
    Y.train <- sc.ref %*% sc.props[, samples$overlapping, 
                                   drop = F]
    X.train <- Biobase::exprs(bulk.eset)[genes, samples$overlapping, 
                                         drop = F]
    X.pred <- Biobase::exprs(bulk.eset)[genes, samples$remaining, 
                                        drop = F]
    template <- base::numeric(base::length(samples$remaining))
    base::names(template) <- samples$remaining
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    Y.pred <- base::matrix(base::vapply(X = genes, FUN = SupervisedTransformBulk, 
                                        FUN.VALUE = template, Y.train, X.train, X.pred, USE.NAMES = TRUE), 
                           nrow = base::length(samples$remaining))
    sample.names <- samples$remaining
  }
  else {
    if (verbose) {
      base::message("Inferring bulk transformation from single-cell alone.")
    }
    Y.train <- sc.ref %*% sc.props
    X.pred <- Biobase::exprs(bulk.eset)[genes, , drop = F]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    Y.pred <- base::matrix(base::vapply(X = genes, FUN = SemisupervisedTransformBulk, 
                                        FUN.VALUE = template, Y.train, X.pred, USE.NAMES = TRUE), 
                           nrow = base::length(sample.names))
  }
  indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
    base::anyNA(column)
  })
  if (base::any(indices)) {
    if (verbose) {
      n.dropped <- base::sum(indices)
      base::message(base::sprintf("Dropped an additional %i genes", 
                                  n.dropped), " for which a transformation could not be learned.")
    }
    if (sum(!indices) == 0) {
      base::stop("Zero genes left for decomposition.")
    }
    Y.pred <- Y.pred[, !indices, drop = F]
    sc.ref <- sc.ref[!indices, , drop = F]
  }
  results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
    sol <- lsei::pnnls(sc.ref, b, sum = 1)
    return(base::append(sol$x, sol$rnorm))
  }))
  base::rownames(results) <- base::append(base::colnames(sc.ref), 
                                          "rnorm")
  base::colnames(results) <- sample.names
  rnorm <- results["rnorm", , drop = T]
  names(rnorm) <- sample.names
  results <- base::list(bulk.props = results[base::colnames(sc.ref), 
                                             , drop = F], sc.props = sc.props, rnorm = rnorm, genes.used = base::rownames(sc.ref))
  return(results)
}

GetOverlappingGenes <- function(sc.eset, bulk.eset, markers, verbose) {
  bulk.genes <- Biobase::featureNames(bulk.eset)
  sc.genes <- Biobase::featureNames(sc.eset)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No overlapping genes found between bulk and ",
                            "single-cell expression."))
  }
  overlapping.genes <- base::intersect(overlapping.genes, markers)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No marker genes found in both bulk and ",
                            "single-cell expression."))
  }
  if (verbose) {
    n.genes <- base::length(overlapping.genes)
    base::message(base::sprintf("Using %i genes in both", n.genes),
                  " bulk and single-cell expression.")
  }
  return(overlapping.genes)
}
CountsToCPM <- function(eset) {
  dropCell <- base::colSums(Biobase::exprs(eset))>200
  message(sum(!dropCell)," cells were removed due to low reads coverage (<200)")
  eset <- eset[,dropCell]
  Biobase::exprs(eset) <- base::sweep(Biobase::exprs(eset),
                                      2, base::colSums(Biobase::exprs(eset)),
                                      `/`) * 1000000
  indices <- base::apply(Biobase::exprs(eset), MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    n.cells <- base::sum(indices)
    base::stop(base::sprintf("Zero expression in selected genes for %i cells",
                             n.cells))
  }
  return(eset)
}
FilterZeroVarianceGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, stats::var) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i zero variance genes.",
                                genes.filtered))
  }
  return(eset)
}
FilterUnexpressedGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, base::sum) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i unexpressed genes.",
                                genes.filtered))
  }
  return(eset)
}
GenerateSCReference <- function(sc.eset, cell.types) {
  cell.labels <- base::factor(sc.eset[[cell.types]])
  all.cell.types <- base::levels(cell.labels)
  aggr.fn <- function(cell.type) {
    base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == cell.type, drop=F])
  }
  template <- base::numeric(base::nrow(sc.eset))
  sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
  return(sc.ref)
}
CalculateSCCellProportions <- function(sc.eset, subject.names, cell.types) {
  individual.labels <- base::factor(sc.eset[[subject.names]])
  individuals <- base::levels(individual.labels)
  cell.labels <- sc.eset[[cell.types]]
  aggr.fn <- function(individual) {
    base::table(cell.labels[individual.labels == individual]) /
      base::length(cell.labels[individual.labels == individual])
  }
  sc.props <- base::sapply(individuals, aggr.fn)
  return(sc.props)
}
SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {
  # Learns linear transformation of observed bulk to match distribution of
  # weighted sum of reference
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- base::scale(Y.train[gene,,drop=T])
  Y.center <- base::attr(Y.train.scaled, "scaled:center")
  Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
  n <- base::length(Y.train.scaled)
  # Shrinkage estimator that minimizes MSE for scaling factor
  shrink.scale <- base::sqrt(base::sum((Y.train[gene,,drop=T]-Y.center)^2)/n+1)
  X.pred.scaled <- base::scale(X.pred[gene,,drop=T])
  Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + Y.center,
                         dimnames=base::list(base::colnames(X.pred), gene))
  return(Y.pred)
}


