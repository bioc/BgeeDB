context("topAnat on propagated and observed calls")

test_that("BgeeDB runs the classic Fisher test of topGO on the propagated calls itself, without propagating them again", {
  data <- smallTopAnatData("propagated")
  object <- topAnat(data, smallGeneList(), nodeSize = 1)
  expect_is(object, "topAnatPropagatedData")
  ## A does not get the genes of its child C
  expect_equal(sort(genesInTerm(object, "A")[["A"]]), sort(paste0("G", 1:10)))
  expect_equal(termStat(object, "A")$Annotated, 10)
  expect_equal(termStat(object, "A")$Significant, 5)
  expect_equal(termStat(object, "A")$Expected, 2.5)
  ## The root holds the 20 genes with expression data, G21 has none
  expect_equal(genesInTerm(object, "BGEE:0")[["BGEE:0"]], sort(paste0("G", 1:20)))
  expect_equal(sigGenes(object), paste0("G", 1:5))
  expect_setequal(graph::nodes(graph(object)), c("A", "B", "C", "R", "BGEE:0"))
  expect_warning(termStat(object, c("A", "D")), "not present")

  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  expect_is(results, "topGOresult")
  ## The one-sided Fisher test of A among the 20 genes with expression data, G21 being left out
  expect_equal(score(results)[["A"]],
               fisher.test(matrix(c(5, 0, 5, 10), nrow = 2), alternative = "greater")$p.value)
  ## Each condition gets the p-value of the contingency table of topGO, 1 without any foreground gene
  for (condition in names(score(results))) {
    annotated <- length(genesInTerm(object, condition)[[1]])
    significant <- length(intersect(genesInTerm(object, condition)[[1]], sigGenes(object)))
    contingency <- cbind(c(significant, 5 - significant), c(annotated - significant, 20 - annotated - 5 + significant))
    expect_identical(score(results)[[condition]],
                     if (significant == 0) 1 else fisher.test(contingency, alternative = "greater")$p.value)
  }
  invisible(capture.output(table <- makeTable(data, object, results, cutoff = 1)))
  expect_equal(table["A", c("annotated", "significant", "expected")], data.frame(annotated = 10L, significant = 5L,
                                                                              expected = 2.5, row.names = "A"))

  ## The node size leaves out the conditions with fewer genes, the root being kept
  expect_setequal(names(score(runTest(topAnat(data, smallGeneList(), nodeSize = 5), "classic", "fisher"))),
                  c("A", "B", "R", "BGEE:0"))
  expect_error(runTest(object, algorithm = "weight", statistic = "fisher"), "observed")
  expect_error(runTest(object, statistic = "fisher"), "observed")
  expect_error(runTest(object, algorithm = "classic", statistic = "ks"), "Fisher")
})

test_that("topGO propagates the observed calls, and only runs its decorrelation algorithms on them from Bgee 16", {
  object <- topAnat(smallTopAnatData("observed"), smallGeneList(), nodeSize = 1)
  ## A gets the genes of its child C
  expect_equal(sort(genesInTerm(object, "A")[["A"]]), sort(paste0("G", c(1:10, 16))))
  expect_error(runTest(object, algorithm = "classic", statistic = "fisher"), "propagated")
  expect_is(runTest(object, algorithm = "weight", statistic = "fisher"), "topGOresult")

  ## The releases prior to 16 only have observed calls, all the algorithms are run on them
  object15 <- topAnat(smallTopAnatData("observed", list(release = "15_2", sendStats = FALSE)),
                      smallGeneList(), nodeSize = 1)
  expect_is(runTest(object15, algorithm = "classic", statistic = "fisher"), "topGOresult")
})
