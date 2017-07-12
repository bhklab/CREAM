
###################### WindowVec is a function to specify window size for each Order of COREs
###################### Input variables of this function are as follows:
###################### 1) InputData: The input data as a table including chromosome regions
###################### in which the first column is chromosome annotation, 
###################### and second and third columns are start and ending positions.
###################### 2) peakNumMin: minimum Order of COREs

WindowVec <- function(InputData, peakNumMin){
  
  WindowVec_Act <- c()
  WindowSize <- WindowSizeRecog(InputData, peakNumMin)
  WindowVec_Act <- c(WindowVec_Act, WindowSize)
  OrderIter <- 1
  while(WindowVec_Act[length(WindowVec_Act)] >0){
    OrderIter <- (OrderIter+1)
    peakNumMax <- OrderIter
    
    WindowSize <- WindowSizeRecog(InputData, peakNumMax)
    
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

