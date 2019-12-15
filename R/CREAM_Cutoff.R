#' CREAM_Cutoff is the function to adapt CREAM for LOCK identification
#'
#' @param in_path Path to the input file (The file inclusing the functional
#' regions)
#' Note. You have to make sure that there is no overlapping regions within the
#' input file
#' @param MinLength Criteria for the minimum number of functional regions in the
#'  input file
#' @param peakNumMin Minimum number of peaks for CORE identification
#' @return Bed file including the identified COREs
#' @examples
#' CREAM(system.file("extdata", "A549_Chr21.bed", package = "CREAM"),
#' MinLength = 1000, peakNumMin = 2)
#' @importFrom utils read.table 
#' @export

CREAM_Cutoff <- function(in_path, MinLength = 1000, peakNumMin = 2){
  
  ClustList <- list()
  MedWidth <- c(1e5)
  SumVec <- c(0)
  MaxOrder <- c(1)
  CutIter <- 1
  for(Cutoff in seq(15,-20,-0.5)){
    CutIter <- (CutIter+1)
    
    aa <- CREAM(in_path,WScutoff = Cutoff*0.1, MinLength = 1000, peakNumMin = 2)
    if(!is.null(aa)){
      ClustList[[CutIter]] <- matrix(aa, ncol = 6)
      if(nrow(ClustList[[CutIter]]) == 0){
        MedWidth <- c(MedWidth, 0)
      }else{
        SumVec <- c(SumVec, sum(as.numeric(ClustList[[CutIter]][,5]))/sum((as.numeric(ClustList[[CutIter]][,3])-
                                                                             as.numeric(ClustList[[CutIter]][,2]))))
      }
    }else{
      ClustList[[CutIter]] <- matrix(rep(NA,3), ncol = 3)
      SumVec <- c(SumVec, 0)
    }
    ##############
    if(length(SumVec) > 2 & min(SumVec[(length(SumVec)-2):(length(SumVec))]) > 0){
      if(SumVec[(length(SumVec)-2)] < 0.95 &
         (3*SumVec[(length(SumVec)-2)]-4*SumVec[(length(SumVec)-1)]
          +SumVec[(length(SumVec))])/SumVec[(length(SumVec)-2)] > 0.05){

        return(ClustList[[(CutIter-2)]])
        stop("Clusters obtained")
      }
    }
  }

  return(ClustList[[1]])

}

