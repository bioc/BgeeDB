context("topAnat on propagated and observed calls")

test_that("topGO does not propagate again the propagated calls, and runs its classic Fisher test on them", {
  object <- topAnat(smallTopAnatData("propagated"), smallGeneList(), nodeSize = 1)
  expect_is(object, "topAnatData")
  ## A does not get the genes of its child C
  expect_equal(sort(genesInTerm(object, "A")[["A"]]), sort(paste0("G", 1:10)))
  expect_equal(termStat(object, "A")$Annotated, 10)
  expect_equal(termStat(object, "A")$Significant, 5)

  results <- runTest(object, algorithm = "classic", statistic = "fisher")
  ## The one-sided Fisher test of A among the 20 genes with expression data, G21 being left out
  expect_equal(score(results)[["A"]],
               fisher.test(matrix(c(5, 0, 5, 10), nrow = 2), alternative = "greater")$p.value)
  expect_error(runTest(object, algorithm = "weight", statistic = "fisher"), "observed")
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
