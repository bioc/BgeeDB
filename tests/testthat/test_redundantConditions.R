context("Conditions redundant with a more precise one")

test_that("a condition is redundant when one of its sub-conditions, at any depth, has the same result", {
  ## G is the parent of I and D, both parents of C; S, sibling of I and D, is not related to C
  relationships <- list(C = c("I", "D"), I = "G", D = "G", S = "G")
  redundantOf <- function(annotated, significant, pValues) {
    ids <- names(annotated)
    setNames(BgeeDB:::.redundantConditions(ids, annotated, significant, pValues[ids], relationships), ids)
  }
  ## The parent with the result of its child is redundant, the child is kept
  expect_equal(redundantOf(c(C = 10, I = 10), c(C = 5, I = 5), c(C = 0.01, I = 0.01)), c(C = FALSE, I = TRUE))
  ## Every ancestor along a chain with the same result, the most precise condition only being kept
  expect_equal(redundantOf(c(C = 10, I = 10, G = 10), c(C = 5, I = 5, G = 5), c(C = 0.01, I = 0.01, G = 0.01)),
               c(C = FALSE, I = TRUE, G = TRUE))
  ## Through an intermediate condition with another result
  expect_equal(redundantOf(c(C = 10, I = 12, G = 10), c(C = 5, I = 5, G = 5), c(C = 0.01, I = 0.02, G = 0.01)),
               c(C = FALSE, I = FALSE, G = TRUE))
  ## The same counts with another p-value, as the decorrelation algorithms can give
  expect_equal(redundantOf(c(C = 10, I = 10), c(C = 5, I = 5), c(C = 0.01, I = 0.5)), c(C = FALSE, I = FALSE))
  ## The same p-value with other counts
  expect_equal(redundantOf(c(C = 10, I = 12), c(C = 5, I = 6), c(C = 0.01, I = 0.01)), c(C = FALSE, I = FALSE))
  ## Two conditions with the same result that are not related
  expect_equal(redundantOf(c(C = 10, S = 10), c(C = 5, S = 5), c(C = 0.01, S = 0.01)), c(C = FALSE, S = FALSE))
})

## The small topAnat data, with D between A and the root R, holding exactly the genes of A
dataWithRedundantParent <- function() {
  data <- smallTopAnatData("propagated")
  data$gene2anatomy <- lapply(data$gene2anatomy, function(conditions) {
    if ("A" %in% conditions) c(conditions, "D") else conditions
  })
  data$organ.relationships$A <- "D"
  data$organ.relationships$D <- "R"
  data$organ.names <- rbind(data$organ.names, data.frame(ID = "D", NAME = "d", stringsAsFactors = FALSE))
  data
}

test_that("makeTable leaves out the redundant conditions by default, and flags them on request", {
  data <- dataWithRedundantParent()
  object <- topAnat(data, smallGeneList(), nodeSize = 1)
  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  invisible(capture.output(all <- makeTable(data, object, results, cutoff = 1, keepRedundant = TRUE)))
  invisible(capture.output(table <- makeTable(data, object, results, cutoff = 1)))

  ## D has the result of A, its child, and the root added by BgeeDB the 20 genes of R, its child
  expect_setequal(all$organId[all$redundant], c("D", "BGEE:0"))
  expect_equal(table$organId, all$organId[!all$redundant])
  expect_false("redundant" %in% names(table))
  ## The FDR is the one of all the conditions tested
  expect_equal(table$FDR, all$FDR[!all$redundant])
  expect_error(makeTable(data, object, results, keepRedundant = NA), "keepRedundant")
})

test_that("printTopAnatGraph draws the most significant conditions that are not redundant", {
  data <- dataWithRedundantParent()
  object <- topAnat(data, smallGeneList(), nodeSize = 1)
  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  ## A is the most significant condition, D repeats its result, C comes next
  expect_identical(score(results)[["D"]], score(results)[["A"]])
  expect_equal(names(BgeeDB:::.mostSignificantConditions(data, object, score(results), 2, FALSE)), c("A", "C"))
  expect_setequal(names(BgeeDB:::.mostSignificantConditions(data, object, score(results), 2, TRUE)), c("A", "D"))

  skip_if_not_installed("Rgraphviz")
  fileName <- tempfile(fileext = ".pdf")
  expect_equal(printTopAnatGraph(data, object, results, firstSigNodes = 2, fileName = fileName), fileName)
  expect_true(file.info(fileName)$size > 0)
  expect_error(printTopAnatGraph(data, object, results, fileName = fileName, keepRedundant = "no"), "keepRedundant")
  unlink(fileName)
})
