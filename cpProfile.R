#!/usr/bin/env Rscript

if(!require(cellmap)){
	stop("Please install the cellmap first!")
}

strPath <- dirname(system.file("extdata","bulk.txt",package="cellmap"))
for(one in list.files("profiles","rds$",full.names=T)){
	file.copy(one,strPath,overwrite=T)
}

message("Successfully copied the profiles. And the examples in the help should work now!")