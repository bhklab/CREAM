#' WindowSizeRecog is a function to specify window size for each order of COREs
#'
#' @param InputData The input data as a table including chromosome regions
#' in which the first column is chromosome annotation, and second and third
#' columns are start and ending positions.
#' @param COREorder Order of the COREs which window size has to be determined for.
#' @param WScutoff Threshold used to identify WS within distribution of maximum distance between peaks for each order of CORE
#' @return Window size identified for each order of CORE
#' @importFrom stats median quantile
#' @examples
#' InputData <- read.table(system.file("extdata", "A549_Chr21.bed",
#' package = "CREAM"), sep="\t")
#' colnames(InputData) <- c("chr", "start", "end")
#' MinLength <- 1000
#' if(nrow(InputData) < MinLength){
#'    stop(paste( "Number of functional regions is less than ", MinLength,
#'    ".", sep = "", collapse = ""))
#' }
#' peakNumMin <- 2
#' WScutoff <- 1.5
#' WindowSize <- WindowSizeRecog(InputData, peakNumMin, WScutoff)
#' @export
WindowSizeRecog <- function(InputData, COREorder, WScutoff){
  
  ChrSeq             <- as.character(unique(InputData[,1]))
  OrderSeqAll_Vec    <- c()
  WindowAll_Vec      <- c()
  
  for(chrIter in ChrSeq){
    
    InputData_Start  <- InputData[which(InputData[,1] == chrIter),"start"]
    InputData_End    <- InputData[which(InputData[,1] == chrIter),"end"]
    RemInd <- which(duplicated(paste(InputData_Start, InputData_End, sep = "_")))
    if(length(RemInd) > 0){
      InputData_Start <- InputData_Start[-RemInd]
      InputData_End <- InputData_End[-RemInd]
    }
    InputData_End    <- InputData_End[order(InputData_Start, decreasing = F)]
    InputData_Start  <- InputData_Start[order(InputData_Start, decreasing = F)]
    InputData_Center <- 0.5*(InputData_Start + InputData_End)
    
    InputData_StartSeq <- min(InputData_Start)
    InputData_EndSeq   <- max(InputData_End)
    
    peakNumIter <- COREorder
    
    if(((length(InputData_Start)-(peakNumIter - 1))-1) > 0){
      OrderElement_Vec  <- rep(0,((length(InputData_Start)-(peakNumIter - 1))-1))
      WindowElement_Vec <- rep(0,((length(InputData_Start)-(peakNumIter - 1))-1))
      ##########
      i <- 1
      while(i < (length(InputData_Start)-(peakNumIter - 1))){
        
        widthElement <- (InputData_End[(i+(peakNumIter - 1))] - InputData_Start[i])
        checkwindow  <- max(InputData_Start[(i+1):(i + (peakNumIter - 1))] - 
                              InputData_End[i:(i+ (peakNumIter - 1) - 1)])
        
        OrderElement_Vec[i]  <- peakNumIter
        WindowElement_Vec[i] <- checkwindow
        i <- i + 1
      }
      
      ##### Window-based analysis
      OrderSeqAll_Vec    <- c(OrderSeqAll_Vec, OrderElement_Vec)
      WindowAll_Vec <- c(WindowAll_Vec, WindowElement_Vec)
    }
  }
  
  i <- COREorder
  SortedWindow_Vec <- sort(WindowAll_Vec[which(OrderSeqAll_Vec == i)])
  
  SortedWindowQuan <- quantile(SortedWindow_Vec)
  aa <- (as.numeric(SortedWindowQuan[4]) + WScutoff*(as.numeric(SortedWindowQuan[4])-
                                                       as.numeric(SortedWindowQuan[2])))
  
  RemovePeaks <- which(SortedWindow_Vec > aa)
  
  if(length(RemovePeaks) > 0){
    bb <- log(SortedWindow_Vec[-RemovePeaks])
  }else{
    bb <- log(SortedWindow_Vec)
  }
  
  bb_quan <- quantile(bb)
  TightReg <- (as.numeric(bb_quan[2]) - WScutoff*(as.numeric(bb_quan[4]) -
                                                    as.numeric(bb_quan[2])))
  Outliers <- which(bb < TightReg)
  
  if(length(Outliers) > 0){
    WindowSize <- exp(TightReg)
  }else{
    WindowSize <- 1
  }
  return(WindowSize)
}
