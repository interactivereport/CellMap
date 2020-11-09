#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(optparse)))

## get input parameters ----------------
opls <- list()
opls[[length(opls)+1]] <- make_option(c('--mixBulk','-b'), action="store", type="character", dest="strBulk", default=NULL,
                                      help="The path to the file (txt:'Tab' separated) which contains mix bulk expression (row:genes; column: samples). [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--feature','-f'), action="store", type="character", dest="strProfile", default=NULL,
                                      help="The full path to the profile file rds for the CellMap method.")
opls[[length(opls)+1]] <- make_option(c('--outPrefix','-o'), action="store", type="character", dest="outPrefix", default=NULL, 
                                      help="The file prefix including path to the output files.")
opls[[length(opls)+1]] <- make_option(c('--delCT','-d'), action="store", type="character", dest="delCT", default=NULL, 
                                      help="a string of removal cell type names separated by ';'.")
opls[[length(opls)+1]] <- make_option(c('--cellCol','-c'), action="store", type="character", dest="cellCol", default=NULL, 
                                      help="a string of cell color names separated by ';'.")
opls[[length(opls)+1]] <- make_option(c('--sigCutoff','-s'), action="store", type="numeric", dest="sigCutoff", default=0.05, 
                                      help="The significance cutoff to call a cell type [default: %default].")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=2,
                                      help="The number of cores to be used [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "exeCellMap.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

if(is.null(args$strBulk) || is.null(args$strProfile)) {
  message("Incomplete arguments list!")
  print_help(opts)
  q()
}
cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"),sep="")
cat("\n")
cat("\n")

source("cellMapDecom.R")
strPrefix <- substr(args$strBulk,1,nchar(args$strBulk)-4)
if(!is.null(args$outPrefix)) strPrefix <- args$outPrefix
delCT <- NULL
if(!is.null(args$delCT) && nchar(args$delCT)>3) delCT <- trimws(unlist(strsplit(args$delCT,";")))
cellCol <- NULL
if(!is.null(args$cellCol) && nchar(args$cellCol)>3){
  tmp <- trimws(unlist(strsplit(args$cellCol,";")))
  cellCol <- setNames(sapply(strsplit(tmp,"="),tail,1),
                      sapply(strsplit(tmp,"="),head,1))
}
print(cellCol)
q()
cellMapDecom(args$strBulk,args$strProfile,strPrefix,delCT,ensemblPath="Data/",pCutoff=args$sigCutoff,cellCol=cellCol)

