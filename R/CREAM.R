#' CREAM is the main function for CORE identification
#'
#' @param in_path Path to the input file (The file inclusing the functional
#' regions)
#' Note. You have to make sure that there is no overlapping regions within the
#' input file
#' @param out_path The path in which you want to store the identified COREs
#' @param MinLength Criteria for the minimum number of functional regions in the
#'  input file
#' @param peakNumMin Minimum number of peaks for CORE identification
#' @param WScutoff Threshold used to identify WS within distribution of maximum
#' distance between peaks for each order of CORE
#' @return Bed file including the identified COREs
#' @examples
#' CREAM(system.file("extdata", "A549_Chr21.bed", package = "CREAM"),
#' system.file("extdata", "A549_Chr21_COREs.bed", package = "CREAM"),
#' MinLength = 1000, peakNumMin = 2)
#' @importFrom utils read.table write.table
#' @export
CREAM <- function(in_path, out_path, WScutoff = 1.5, MinLength = 1000, peakNumMin = 2){
  InputData   <- read.table(in_path, sep="\t")
  colnames(InputData) <- c("chr", "start", "end")
  ###########################
  ###########################
  if(nrow(InputData) < MinLength){
    stop(paste( "Number of functional regions is less than ", MinLength, ".", sep = "", collapse = ""))
  }
  #####################
  WindowVecFinal <- WindowVec(InputData, peakNumMin, WScutoff)
  OutputList <- ElementRecog(InputData, WindowVecFinal, (1+length(WindowVecFinal)), peakNumMin)
  WidthSeq_Vec    <-  OutputList[[1]]
  StartSeq_Vec    <-  OutputList[[2]]
  EndSeq_Vec      <-  OutputList[[3]]
  ChrSeq_Vec      <-  OutputList[[4]]
  OrderSeq_Vec    <-  OutputList[[5]]
  WinSizeSeq_Vec  <-  OutputList[[6]]
  ####################
  if(!is.null(StartSeq_Vec)){
    CombinedData <- cbind(ChrSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]],
                          StartSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]],
                          EndSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]],
                          OrderSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]],
                          WidthSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]],
                          WinSizeSeq_Vec[sort(OrderSeq_Vec, decreasing = T, index.return = T)[[2]]])
    colnames(CombinedData) <- c("Chr", "Start", "End", "Order", "Width", "WindowSize")

    MinPeaks <- PeakMinFilt(CombinedData, WindowVecFinal)

    RemovePeaks <- which(as.numeric(CombinedData[,"Order"]) < MinPeaks)
    if(length(RemovePeaks) > 0){
      CombinedData <- CombinedData[-RemovePeaks,]
    }

    CombinedData <- CombinedData[,c(1:3)]
    colnames(CombinedData) <- NULL
    write.table(CombinedData , file = out_path,
                row.names=FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
}
