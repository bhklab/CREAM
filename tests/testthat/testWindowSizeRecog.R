library(CREAM)

context("WindowSizeRecog")

test_that("Specify window size for each order of COREs", {
    InputData   <- read.table("A549_Chr21.bed", sep="\t")
    colnames(InputData) <- c("chr", "start", "end")
    MinLength <- 1000
    if(nrow(InputData) < MinLength){
        stop(paste( "Number of functional regions is less than ", MinLength, ".", sep = "", collapse = ""))
    }
    peakNumMin <- 2
    WScutoff <- 1.5
    WindowSize <- WindowSizeRecog(InputData, peakNumMin, WScutoff)
    expect_equal_to_reference(WindowSize, "windowSize.rds")
})