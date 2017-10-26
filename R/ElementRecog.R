#' ElementRecog is a function to identify COREs
#'
#' @param InputData The input data as a table including chromosome regions
#' in which the first column is chromosome annotation, and second and third
#' columns are start and ending positions.
#' @param windowSize_Vec Vector of window sizes ordered based on order of CORE
#' @param peakNumMax Maximum order of COREs (e.g. maximum number of peaks within COREs)
#' @param peakNumMin Minimum order of COREs (e.g. minimum number of peaks within COREs)
#' @return Identified COREs for the given input regions
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
#' WindowVecFinal <- WindowVec(InputData, peakNumMin, WScutoff)
#' OutputList <- ElementRecog(InputData, WindowVecFinal,
#' (1+length(WindowVecFinal)), peakNumMin)
#' @export
ElementRecog <- function(InputData, windowSize_Vec, peakNumMax, peakNumMin){

  ChrSeq             <- as.character(unique(InputData[,1]))
  WidthSeq_All       <- c()
  StartRegionAll_Vec <- c()
  EndRegionAll_Vec   <- c()
  ChrSeqAll_Vec      <- c()
  OrderSeqAll_Vec    <- c()
  SDSeqAll_Vec       <- c()
  WindowSizeAll_Vec  <- c()

  for(chrIter in ChrSeq){

    InputData_Start  <- InputData[which(InputData[,1] == chrIter),"start"]
    InputData_End    <- InputData[which(InputData[,1] == chrIter),"end"]

    InputData_End    <- InputData_End[order(InputData_Start, decreasing = F)]
    InputData_Start  <- InputData_Start[order(InputData_Start, decreasing = F)]
    InputData_Center <- 0.5*(InputData_Start + InputData_End)

    InputData_StartSeq <- min(InputData_Start)
    InputData_EndSeq   <- max(InputData_End)

    ChrElement_Vec    <- c()
    StartElement_Vec  <- c()
    EndElement_Vec    <- c()
    WidthElement_Vec  <- c()
    OrderElement_Vec  <- c()
    SDElement_Vec     <- c()
    WindowSize_Vec <- c()

    for(peakNumIter in seq(peakNumMax, peakNumMin, by = -1)){
      i <- 1
      WindowSize <- windowSize_Vec[(peakNumIter - 1)]
      while(i < (length(InputData_Start)-(peakNumIter - 1))){

        widthElement <- (InputData_End[(i+(peakNumIter - 1))] - InputData_Start[i])
        checkwindow  <- max(InputData_Start[(i+1):(i + (peakNumIter - 1))] -
                              InputData_End[i:(i+ (peakNumIter - 1) - 1)])

        if(checkwindow < WindowSize){
          ChrElement_Vec    <- c(ChrElement_Vec, chrIter)
          StartElement_Vec  <- c(StartElement_Vec, InputData_Start[i])
          EndElement_Vec    <- c(EndElement_Vec, InputData_End[(i+(peakNumIter-1))])
          WidthElement_Vec  <- c(WidthElement_Vec, widthElement)
          OrderElement_Vec  <- c(OrderElement_Vec, peakNumIter)
          WindowSize_Vec <- c(WindowSize_Vec, checkwindow)
          InputData_Start <- InputData_Start[-(i:(i+(peakNumIter-1)))]
          InputData_End   <- InputData_End[-(i:(i+(peakNumIter-1)))]
          InputData_Center <- InputData_Center[-(i:(i+(peakNumIter-1)))]
        }else{
          i <- i + 1
        }
      }
    }

    ##### Window-based analysis
    WidthSeq_All       <- c(WidthSeq_All, WidthElement_Vec)
    StartRegionAll_Vec <- c(StartRegionAll_Vec, StartElement_Vec)
    EndRegionAll_Vec   <- c(EndRegionAll_Vec, EndElement_Vec)
    ChrSeqAll_Vec      <- c(ChrSeqAll_Vec, ChrElement_Vec)
    OrderSeqAll_Vec    <- c(OrderSeqAll_Vec, OrderElement_Vec)
    WindowSizeAll_Vec  <- c(WindowSizeAll_Vec, WindowSize_Vec)
  }

  return(list(WidthSeq_All,  StartRegionAll_Vec,
              EndRegionAll_Vec, ChrSeqAll_Vec, OrderSeqAll_Vec,
              WindowSizeAll_Vec))

}
