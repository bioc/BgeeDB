context("cell processed data")

## These tests do not download any H5AD file: the arguments are checked before the download

test_that("getCellProcessedData rejects an unsupported R package before downloading anything", {
  bgee <- Bgee$new(species = "Gallus_gallus", dataType = "sc_droplet_based")
  expect_error(getCellProcessedData(bgee, experimentId = "ERP132576", package = "foo"),
               "should either be anndata or zellkonverter")
  expect_error(getCellProcessedData(bgee, experimentId = "ERP132576",
                                    package = c("anndata", "zellkonverter")),
               "should either be anndata or zellkonverter")
})

test_that("Galaxy wrappers can build a Bgee object without datatype, and downloading cell data needs a single cell datatype", {
  ## the datatype is not mandatory to build the object used by the wrappers
  expect_error(createBgeeObject_galaxy(species = "Gallus_gallus"), NA)
  ## but a single cell datatype is needed to download H5AD files: without any datatype the object
  ## holds all of them, which has to give the explicit message and not an error on the condition
  expect_error(downloadCellProcessedData_galaxy(species = "Gallus_gallus", experimentId = "ERP132576",
                                                cellProcessedDataFile = tempfile()),
               "only available for single cell data")
  expect_error(downloadCellProcessedData_galaxy(species = "Gallus_gallus", dataType = "rna_seq",
                                                experimentId = "ERP132576",
                                                cellProcessedDataFile = tempfile()),
               "only available for single cell data")
})
