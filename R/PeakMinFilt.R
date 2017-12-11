
#' PeakMinFilt is a function to filter the lowest Order of COREs which distance
#' between functional regions is close to the corresponding Window Size
#'
#' @param Clusters_init Table of indetified COREs before filteration
#' @param WindowVecFinal Vector of window sizes ordered based on order of CORE
#' @return Minimum order of COREs
#' @importFrom stats median
#' @export
PeakMinFilt <- function(Clusters_init, WindowVecFinal){
  print("Filtering clusters of low order")
  UniqueOrder <- sort(unique(as.numeric(Clusters_init[,"Order"])), decreasing = F)
  Zscore <- c()
  for(OrderIter in 1:length(UniqueOrder)){
    TargetInd <- which(as.numeric(Clusters_init[,"Order"]) == UniqueOrder[OrderIter])
    Zscore <- c(Zscore, (WindowVecFinal[(UniqueOrder[OrderIter] - 1)] - 
                           median(as.numeric(Clusters_init[TargetInd,"WindowSize"]))
    )/median(as.numeric(Clusters_init[TargetInd,"WindowSize"])))
    
  }
  
  TargetOrder <- 0
  if(length(Zscore) > 2){
    for(OrderIter in 2:length(Zscore)){
      if(Zscore[OrderIter] < Zscore[(OrderIter-1)]){
        return(UniqueOrder[OrderIter])
        stop("Vector obtained")
      }
    }
  }
  return(max(UniqueOrder))
}
