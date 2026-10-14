context("Data types of a release")

test_that("A Bgee object created without any data type does not use Affymetrix and EST data from Bgee 16", {
  expect_equal(.releaseDataTypes("15_2"),
               c("rna_seq", "affymetrix", "est", "in_situ", "sc_full_length", "sc_droplet_based"))
  expect_equal(.releaseDataTypes("16_0"), c("rna_seq", "in_situ", "sc_full_length", "sc_droplet_based"))
  ## The release can be written with a dot
  expect_equal(.releaseDataTypes("16.0"), .releaseDataTypes("16_0"))
})
