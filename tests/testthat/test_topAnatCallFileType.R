context("topAnat on propagated and observed calls")

## 20 genes expressed at the root R; A: genes 1 to 10; B: genes 6 to 15; C, child of A: genes 1 to 3
## and gene 16, absent from A: the propagated calls of a condition do not always hold those of its children
smallTopAnatData <- function(callFileType, bgeeObject = NULL) {
  genes <- paste0("G", 1:20)
  conditions <- list(A = genes[1:10], B = genes[6:15], C = genes[c(1:3, 16)])
  gene2anatomy <- lapply(genes, function(g) c("R", names(conditions)[sapply(conditions, function(c) g %in% c)]))
  names(gene2anatomy) <- genes
  list(gene2anatomy = gene2anatomy,
       organ.relationships = list(A = "R", B = "R", C = "A", R = "BGEE:0"),
       organ.names = data.frame(ID = c("R", "A", "B", "C", "BGEE:0"),
                                NAME = c("root", "a", "b", "c", "root"), stringsAsFactors = FALSE),
       bgee.object = bgeeObject, call.file.type = callFileType)
}
## Foreground: genes 1 to 5, and G21 that has no expression data
smallGeneList <- function() {
  geneList <- factor(c(rep(1, 5), rep(0, 15), 1))
  names(geneList) <- paste0("G", 1:21)
  geneList
}

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
