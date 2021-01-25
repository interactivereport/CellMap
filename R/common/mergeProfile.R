###############################
## mergeProfile.R
##
###############################

mergeProfile <- function(preProfile,curProfile){
  if(is.null(preProfile)){
    return(curProfile)
  }else{
    Exp <- merge(preProfile$expr,curProfile$expr,all=T,by="row.names",sort=F)
    rownames(Exp) <- Exp[,1]
    Exp <- Exp[,-1]
    Exp[is.na(Exp)] <- 0
    cat(ncol(preProfile$sets),"\t",ncol(curProfile$sets),"\n")
    selG <- list()
    for(i in names(preProfile$selG)) selG[[i]] <- unique(c(preProfile$selG[[i]],curProfile$selG[[i]]))
    return(list(expr=Exp,selG=selG,sets=rbind(preProfile$sets,curProfile$sets)))
  }
}