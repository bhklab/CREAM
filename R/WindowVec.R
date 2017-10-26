#' WindowVec is a function to specify window size for each order of COREs
#'
#' @param InputData The input data as a table including chromosome regions in
#' which the first column is chromosome annotation,  and second and third
#' columns are start and ending positions.
#' @param peakNumMin Minimum order of COREs
#' @param WScutoff Threshold used to identify WS within distribution of maximum distance between peaks for each order of CORE
#' @return Vector of window sizes from order 2 up to maximum order of COREs
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
#' @export
WindowVec <- function(InputData, peakNumMin, WScutoff){

  WindowVec_Act <- c()
  WindowSize <- WindowSizeRecog(InputData, peakNumMin, WScutoff)
  WindowVec_Act <- c(WindowVec_Act, WindowSize)
  OrderIter <- 2
  while(WindowVec_Act[length(WindowVec_Act)] >0){
    OrderIter <- (OrderIter+1)
    peakNumMax <- OrderIter

    WindowSize <- WindowSizeRecog(InputData, peakNumMax, WScutoff)

    WindowVec_Act <- c(WindowVec_Act, WindowSize)
    lenWinVec <- length(WindowVec_Act)

    if(WindowVec_Act[lenWinVec] > 1 & WindowVec_Act[(lenWinVec - 1)] > 1){
      if(WindowVec_Act[lenWinVec] < WindowVec_Act[(lenWinVec - 1)] ||
         WindowVec_Act[lenWinVec] == WindowVec_Act[(lenWinVec - 1)]){
        return(WindowVec_Act[1:(lenWinVec - 1)])
        stop("Vector obtained")
      }
    }
  }
  return(WindowVec_Act)
}

