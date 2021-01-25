###############################
## updateProfile.R
##
##################

updateProfile <- function(allProfile,evalV,para){
  cat("removing",length(evalV$rmSetsIndex),"sets from total of",nrow(allProfile$sets),"sets\n")
  
  if(length(evalV$rmSetsIndex)>0)
    allProfile$sets <- allProfile$sets[-sample(evalV$rmSetsIndex,
                                               min(length(evalV$rmSetsIndex),
                                                   round(para$maxRMrate*nrow(allProfile$sets)))),]

  allProfile$expr <- allProfile$expr[,colnames(allProfile$expr)%in%as.vector(allProfile$sets)]
  para$cellMap <- para$trainMap[para$trainMap%in%evalV$poorCellType]
  allProfile$para <- para
  allProfile$score <- evalV$score
  return(allProfile)
}