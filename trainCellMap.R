#!/usr/bin/env Rscript
suppressMessages(require(optparse,quietly=T, warn.conflicts=F))
## input ----
opls <- list()
opls[[length(opls)+1]] <- make_option(c('--Data','-d'), action="store", type="character", dest="data", default=NULL,
                                      help="The full path to the file listed training sc/sn RNAseq data (required). Three columns separated by tab without header.\n
                                      \tfirst column: the name of the dataset;\n\tsecond column: the full path to the data rds file;\n\t
                                      third column:the cell types defined in the datasets (first element separated by | in column names), e.g. A, B, C,
                                      and the cell type names prepered, e.g. CT1, CT2. In format of: A=CT1;B=CT2;C=CT2.")
opls[[length(opls)+1]] <- make_option(c('--DEG','-g'), action="store", type="character", dest="selFeature", default='edgeR',
                                      help="The method (DESeq2 or edgeR) to extract the cell type genes. [default: %default].")
opls[[length(opls)+1]] <- make_option(c('--bulk','-b'), action="store", type="character", dest="strBulk",default="",
                                      help="The full path to the real bulk expression matrix with known composition (tab separated, row is the genes, column is the samples)")
opls[[length(opls)+1]] <- make_option(c('--rate','-r'), action="store", type="character", dest="strRate",default="",
                                      help="The full path to the cell type compostion matrix (tab separated, row is the cell type, column is the samples).")
opls[[length(opls)+1]] <- make_option(c('--out','-o'), action="store", type="character", dest="output",default=NULL,
                                      help="The output CellMap feature files, including the prefix.")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=8, help="Path to write output files to [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "trainCellMap.R", description = "", epilogue = "")

para <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)
if(is.null(para$data) || is.null(para$output)) {
  message("Incomplete arguments list!")
  print_help(opts)
  q()
}
message(paste(paste(names(para),para,sep=": "),collapse="\n"))

## call the training process ----------
source("cellMapTraining.R")
D <- read.table(para$data,sep="\t",as.is=T,row.names=1)
strD <- setNames(D[,1],rownames(D))
## get the celltypeMap
tmp <- unlist(strsplit(D[,2],";"))
cellTypeID <- sapply(strsplit(tmp,"="),head,1)
cellTypeName <- setNames(sapply(strsplit(tmp,"="),tail,1)[!duplicated(cellTypeID)],cellTypeID[!duplicated(cellTypeID)])

cellMapTraining(strD,paste0(para$output,"/Profile"),cellTypeMap=cellTypeName,core=para$core)


