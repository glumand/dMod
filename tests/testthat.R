library(testthat)
library(dMod)

# Several test files setwd(tempdir()) without restoring (test-petab.R,
# test-deriv2-*, test-noparam-trafos.R, test-reparam.R, test-steadystates.R).
# That breaks testthat's relative-path lookups in later tests. Snapshot the
# original wd here and restore it on exit so test_check's path resolution
# stays consistent across the whole run.
.tt_initial_wd <- getwd()
on.exit(setwd(.tt_initial_wd), add = TRUE)

test_check("dMod")
