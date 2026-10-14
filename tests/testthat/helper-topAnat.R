## Small topAnat data shared by the tests of the topAnat functions

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
