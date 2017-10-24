library(CREAM)

context("CREAM")

test_that("COREs are identifiable", {
    list.files( path="." )
    CREAM("../../inst/testdata/A549_Chr21.bed",
          "../../inst/testdata/out_COREs.bed",
          MinLength = 1000, peakNumMin = 2)
    expect_true("../../inst/testdata/A549_Chr21_COREs.bed",
                "../../inst/testdata/out_COREs.bed")
})