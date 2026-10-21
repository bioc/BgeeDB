context("gene list and background of topAnat")

test_that("the gene list is restricted to a background given as gene IDs", {
  object <- topAnat(smallTopAnatData("propagated"), smallGeneList(), nodeSize = 1, background = paste0("G", 1:10))
  ## G11 to G21 are left out: 10 genes with expression data remain, with the 5 foreground genes
  expect_equal(sum(object@feasible), 10)
  expect_equal(sigGenes(object), paste0("G", 1:5))
  expect_error(topAnat(smallTopAnatData("propagated"), smallGeneList(), nodeSize = 1, background = paste0("G", 6:20)),
               "none of the foreground genes")
})

test_that("from Bgee 16, the gene list is restricted by default to the protein coding genes of Bgee", {
  local_mocked_bindings(.loadTopAnatBackground = function(myBgeeObject, background) {
    expect_equal(background, "protein_coding")
    paste0("G", 1:15)
  }, .package = "BgeeDB")
  bgee16 <- list(release = "16_0", sendStats = FALSE)
  object <- topAnat(smallTopAnatData("propagated", bgee16), smallGeneList(), nodeSize = 1)
  expect_equal(sum(object@feasible), 15)
  object <- topAnat(smallTopAnatData("propagated", bgee16), foregroundGenes = paste0("G", 1:5), nodeSize = 1)
  expect_equal(sum(object@feasible), 15)
  ## A background downloaded from Bgee needs a Bgee object, and a release from Bgee 16
  expect_error(topAnat(smallTopAnatData("propagated"), smallGeneList(), nodeSize = 1, background = "polyA"), "Bgee object")
  expect_error(topAnat(smallTopAnatData("observed", list(release = "15_2", sendStats = FALSE)), smallGeneList(),
                       nodeSize = 1, background = "polyA"), "Bgee 16")
})

test_that("the gene list is not restricted by default before Bgee 16, nor without a Bgee object", {
  expect_equal(sum(topAnat(smallTopAnatData("propagated"), smallGeneList(), nodeSize = 1)@feasible), 20)
  object15 <- topAnat(smallTopAnatData("observed", list(release = "15_2", sendStats = FALSE)), smallGeneList(), nodeSize = 1)
  expect_equal(sum(object15@feasible), 20)
  ## The genes of interest alone are then compared to all the genes with expression data
  object15 <- topAnat(smallTopAnatData("observed", list(release = "15_2", sendStats = FALSE)),
                      foregroundGenes = paste0("G", 1:5), nodeSize = 1)
  expect_equal(sum(object15@feasible), 20)
  expect_equal(sort(sigGenes(object15)), paste0("G", 1:5))
})

test_that("the gene list is built from the genes of interest and the background", {
  data <- smallTopAnatData("propagated")
  ## A custom background: its genes with expression data, G21 having none
  object <- topAnat(data, foregroundGenes = paste0("G", 1:5), background = c(paste0("G", 1:12), "G21"), nodeSize = 1)
  expect_equal(sum(object@feasible), 12)
  expect_equal(sigGenes(object), paste0("G", 1:5))
  ## The same analysis as with the gene list built by hand
  geneList <- factor(as.integer(paste0("G", 1:12) %in% paste0("G", 1:5)))
  names(geneList) <- paste0("G", 1:12)
  expect_identical(score(runTest(object, "classic", "fisher")),
                   score(runTest(topAnat(data, geneList, nodeSize = 1), "classic", "fisher")))
  ## Without a Bgee object nor a background: all the genes with expression data, G21 being left out
  object <- topAnat(data, foregroundGenes = c(paste0("G", 1:5), "G21"), nodeSize = 1)
  expect_equal(sum(object@feasible), 20)
  expect_equal(sigGenes(object), paste0("G", 1:5))
  expect_error(topAnat(data, foregroundGenes = "G21", nodeSize = 1), "none of the genes of interest")
})

test_that("either a gene list or the genes of interest is required, and a gene list must be named", {
  data <- smallTopAnatData("propagated")
  expect_error(topAnat(data, nodeSize = 1), "either a gene list")
  expect_error(topAnat(data, smallGeneList(), foregroundGenes = "G1", nodeSize = 1), "either a gene list")
  expect_error(topAnat(data, c("G1", "G2"), nodeSize = 1), "right format")
})
