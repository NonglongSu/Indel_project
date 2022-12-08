#setwd("~/Dropbox (ASU)/Indel_project/Script/tests")
suppressPackageStartupMessages(library(covr))


covr = file_coverage(source_files = c("../sw_gap.R", 
                                      "../update_gap.R"),
                     test_files   = c("testthat/test_sw_gap.R",
                                      "testthat/test_update_gap.R"))
covr
report(covr,file="coverage_report.html")


#testthat::test_dir("testthat")
#Error: invalid version specification ‘0.68’






################################################################
# cat("add <- function(x, y) { x + y }", file="add.R")
# cat("add(1, 2) == 3", file="add_test.R")
# file_coverage(source_files = "add.R", test_files = "add_test.R")
# file.remove(c("add.R", "add_test.R"))
