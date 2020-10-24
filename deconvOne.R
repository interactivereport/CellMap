#!/usr/bin/env Rscript
## deconvOne.R

suppressMessages(suppressWarnings(require(optparse)))
## get input parameters ----------------
opls <- list()

opls[[length(opls)+1]] <- make_option(c('--mixBulk','-b'), action="store", type="character", dest="strBulk", default=NULL, 
                                      help="The path to the file (txt:'Tab' separated) which contains mix bulk expression (row:genes; column: samples). [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--sc','-f'), action="store", type="character", dest="strSC", default='GSE103723',
                                      help="accession ID or path to the sc/sn RNA expression rds data. The sc/sn RNAseq accession ID,available: GSE76381, GSE103723, phs001836 and GSE67835; the sc/sn RNA expression data should be rds matrix, with genes as row, cell as columns (column names separated by '|' with the first element to be cell type) [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--cutoff','-c'), action="store", type="numeric", dest="cutoff", default=1e6, 
                                      help="Numeric: the total reads/measurements of all cells for a celltype to be considered [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--method','-m'), action="store", type="character", dest="method", default="CellMap", 
                                      help="The methods to be used to decompose (MuSiC or CellMap) [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--outPrefix','-o'), action="store", type="character", dest="prefix", default=paste(getwd(),"Decomposition_",sep=""), 
                                      help="prefix of the result file name [default: %default]")
opls[[length(opls)+1]] <- make_option(c('--thread','-t'), action="store", type="numeric", dest="core", default=8, help="Numeric: the number of cores to be used [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "deconvOne.R", description = "", epilogue = "")

para <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

if(is.null(para$strBulk) || is.null(para$strSC)){
  message("Incomplete arguments list!")
  print_help(opts)
  q()
}
### ----
source("cellMapOne.R")
cellMapOne(para$strBulk,para$strSC,para$prefix,para$method,expCutoff=para$cutoff,core=para$core)


