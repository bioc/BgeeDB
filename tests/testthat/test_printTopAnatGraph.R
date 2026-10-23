context("printTopAnatGraph")

test_that("the graph of the significant conditions is written in a PDF file", {
  skip_if_not_installed("Rgraphviz")
  data <- smallTopAnatData("propagated")
  object <- topAnat(data, smallGeneList(), nodeSize = 1)
  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  fileName <- tempfile(fileext = ".pdf")
  rgraphvizAttached <- "package:Rgraphviz" %in% search()

  expect_equal(printTopAnatGraph(data, object, results, firstSigNodes = 2, fileName = fileName), fileName)
  expect_true(file.info(fileName)$size > 0)
  ## Rgraphviz is only attached while the graph is drawn
  expect_equal("package:Rgraphviz" %in% search(), rgraphvizAttached)
  unlink(fileName)
})

test_that("the graph is drawn when a condition has a p-value of 0", {
  skip_if_not_installed("Rgraphviz")
  data <- smallTopAnatData("propagated")
  object <- topAnat(data, smallGeneList(), nodeSize = 1)
  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  ## The strongest enrichments give a p-value of 0 to the most significant conditions
  scores <- score(results)
  scores[which.min(scores)] <- 0
  score(results) <- scores
  fileName <- tempfile(fileext = ".pdf")

  expect_equal(printTopAnatGraph(data, object, results, firstSigNodes = 2, fileName = fileName), fileName)
  expect_true(file.info(fileName)$size > 0)
  unlink(fileName)
})
