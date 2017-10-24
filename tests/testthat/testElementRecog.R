library(CREAM)

context("ElementRecog")

test_that("Identify COREs", {
    InputData   <- read.table("A549_Chr21.bed", sep="\t")
    colnames(InputData) <- c("chr", "start", "end")
    MinLength <- 1000
    if(nrow(InputData) < MinLength){
        stop(paste( "Number of functional regions is less than ", MinLength, ".", sep = "", collapse = ""))
    }
    peakNumMin <- 2
    WScutoff <- 1.5
    WindowVecFinal <- WindowVec(InputData, peakNumMin, WScutoff)
    OutputList <- ElementRecog(InputData, WindowVecFinal, (1+length(WindowVecFinal)), peakNumMin)
    expect_equal_to_reference(OutputList, "outputList.rds")
})