## The methods of the topAnatPropagatedData class, defined in AllClasses.R: the classic Fisher test on the expression
## calls propagated by Bgee, without topGO.

######################## topAnatPropagatedData class, defined in AllClasses.R ########################
## Builds the topAnatPropagatedData object. The genes and the conditions are selected as topGO selects them:
## the root of the conditions holds every gene with expression data, and the conditions with fewer genes than
## the node size are left out, unless it is 1.
.topAnatPropagatedData <- function(description, ontology, geneList, nodeSize, parentMapping, gene2Nodes) {
  allGenes <- names(geneList)
  ## The selection of topGO for a gene list given as a factor
  selected <- as.logical(as.integer(levels(geneList)))[geneList]
  parentMapping <- as.list(parentMapping)
  nodeSize <- as.integer(max(nodeSize, 1))

  cat("\nCounting the genes of each condition.... ")
  ## The genes of the gene list with expression data, as .annFUN.gene2Nodes selects them for topGO. The list
  ## is an array when loadTopAnatData builds it with tapply: only its names are kept
  genes <- intersect(names(gene2Nodes), allGenes)
  gene2Nodes <- gene2Nodes[genes]
  attributes(gene2Nodes) <- list(names = genes)
  if (any(is.na(gene2Nodes)))
    gene2Nodes <- gene2Nodes[!is.na(gene2Nodes)]
  gene2Nodes <- gene2Nodes[lengths(gene2Nodes) > 0]
  feasible <- allGenes %in% names(gene2Nodes)
  foregroundGenes <- allGenes[selected & feasible]

  ## The numbers of genes and of foreground genes of each condition, in one pass over the gene-condition pairs
  conditions <- unlist(gene2Nodes, use.names = FALSE)
  foreground <- rep(names(gene2Nodes) %in% foregroundGenes, lengths(gene2Nodes))
  conditionIds <- unique(conditions)
  index <- match(conditions, conditionIds)
  annotated <- tabulate(index, length(conditionIds))
  significant <- tabulate(index[foreground], length(conditionIds))
  rm(conditions, foreground, index)

  ## The conditions of the graph: those with genes and all their ancestors, the root holding every gene
  root <- .rootCondition(parentMapping)
  if (length(root) != 1) {
    stop("ERROR: the conditions of the topAnatData object do not have a single root.")
  }
  nodes <- sort(.conditionsAndAncestors(conditionIds, parentMapping))
  nodeIndex <- match(nodes, conditionIds)
  nodeAnnotated <- annotated[nodeIndex]
  nodeAnnotated[is.na(nodeIndex)] <- 0L
  nodeSignificant <- significant[nodeIndex]
  nodeSignificant[is.na(nodeIndex)] <- 0L
  nodeAnnotated[nodes == root] <- length(gene2Nodes)
  nodeSignificant[nodes == root] <- sum(names(gene2Nodes) %in% foregroundGenes)
  counts <- data.frame(Annotated = nodeAnnotated, Significant = nodeSignificant, row.names = nodes)
  ## No condition is left out for a node size of 1, as with topGO
  if (nodeSize > 1) {
    counts <- counts[counts$Annotated >= nodeSize, , drop = FALSE]
  }
  cat(" ( ", sum(feasible), " genes with expression data, ", nrow(counts), " conditions tested. )\n")

  new("topAnatPropagatedData", description = description, ontology = ontology, allGenes = allGenes,
      selected = selected, feasible = feasible, nodeSize = nodeSize, gene2Nodes = gene2Nodes,
      parentMapping = parentMapping, counts = counts)
}

## The root of the conditions: the parent that is the child of no condition
.rootCondition <- function(parentMapping) {
  setdiff(unique(unlist(parentMapping, use.names = FALSE)), names(parentMapping))
}

## The conditions and all their ancestors, up to the root
.conditionsAndAncestors <- function(conditions, parentMapping) {
  nodes <- unique(conditions)
  newNodes <- nodes
  while (length(newNodes) > 0) {
    newNodes <- setdiff(unique(unlist(parentMapping[newNodes], use.names = FALSE)), nodes)
    nodes <- c(nodes, newNodes)
  }
  nodes
}

## The requested conditions that were tested, the others being left out with the warning of topGO for the
## nodes absent from its graph
.testedConditions <- function(object, whichGO) {
  tested <- whichGO %in% rownames(object@counts)
  if (!all(tested)) {
    warning("Nodes not present in the graph:", whichGO[!tested])
    whichGO <- whichGO[tested]
  }
  whichGO
}

## The classic Fisher test of topGO, computed for all the conditions at once from their numbers of genes: same
## contingency table as topGO, and same p-value as fisher.test(alternative = "greater"), that is
## phyper(x - 1, m, n, k, lower.tail = FALSE). As with topGO, the conditions without foreground gene get 1
## without being tested, and come after the others.
#' @rdname topAnat
#' @export
setMethod("runTest", signature(object = "topAnatPropagatedData", algorithm = "character", statistic = "character"),
          function(object, algorithm, statistic, ...) {
            if (tolower(algorithm) != "classic") {
              stop(paste0("ERROR: the expression calls are propagated by Bgee, only the classic algorithm, ",
                          "without decorrelation, can be run on them. The decorrelation algorithms need the observed ",
                          "calls: please retrieve them with loadTopAnatData(..., callFileType = \"observed\")."))
            }
            if (tolower(statistic) != "fisher") {
              stop("ERROR: only the Fisher test can be run on the expression calls propagated by Bgee.")
            }
            message("Running the classic Fisher test on the expression calls propagated by Bgee, without topGO, ",
                    "that would propagate them again.")
            counts <- object@counts
            feasibleGenes <- object@allGenes[object@feasible]
            ## The foreground genes among the genes with expression data, counted as topGO counts them
            numGenes <- length(feasibleGenes)
            numSigAll <- sum(feasibleGenes %in% sigGenes(object))
            nontrivial <- counts$Significant > 0
            pValues <- phyper(counts$Significant[nontrivial] - 1, numSigAll, numGenes - numSigAll,
                              counts$Annotated[nontrivial], lower.tail = FALSE)
            scores <- c(setNames(pValues, rownames(counts)[nontrivial]),
                        setNames(rep(1, sum(!nontrivial)), rownames(counts)[!nontrivial]))
            testedAlgorithm <- "classic"
            attr(testedAlgorithm, "testClass") <- "classicCount"
            new("topGOresult", description = paste(object@description, "\nOntology:", object@ontology, sep = " "),
                score = scores, testName = "fisher", algorithm = testedAlgorithm,
                geneData = c(Annotated = numGenes, Significant = length(sigGenes(object)),
                             NodeSize = object@nodeSize, SigTerms = sum(nontrivial)))
          })

## Without an algorithm, topGO runs weight01, refused on the propagated calls
#' @rdname topAnat
#' @export
setMethod("runTest", signature(object = "topAnatPropagatedData", algorithm = "missing", statistic = "character"),
          function(object, algorithm, statistic, ...) {
            runTest(object, "weight01", statistic, ...)
          })

## The numbers of genes, of foreground genes, and the expected number of foreground genes of conditions, as
## termStat of topGO
#' @rdname topAnat
#' @export
setMethod("termStat", signature(object = "topAnatPropagatedData", whichGO = "character"),
          function(object, whichGO) {
            if (length(whichGO) == 0) {
              return(data.frame(Annotated = integer(), Significant = integer(), Expected = numeric()))
            }
            whichGO <- .testedConditions(object, whichGO)
            counts <- object@counts[whichGO, , drop = FALSE]
            expected <- length(sigGenes(object)) / sum(object@feasible)
            data.frame(Annotated = counts$Annotated, Significant = counts$Significant,
                       Expected = round(counts$Annotated * expected, 2), row.names = whichGO)
          })

#' @rdname topAnat
#' @export
setMethod("termStat", signature(object = "topAnatPropagatedData", whichGO = "missing"),
          function(object, whichGO) {
            termStat(object, rownames(object@counts))
          })

## The genes of conditions, sorted as topGO lists them, the root holding every gene with expression data
#' @rdname topAnat
#' @export
setMethod("genesInTerm", signature(object = "topAnatPropagatedData", whichGO = "character"),
          function(object, whichGO) {
            whichGO <- .testedConditions(object, whichGO)
            conditions <- unlist(object@gene2Nodes, use.names = FALSE)
            inConditions <- conditions %in% whichGO
            geneIndex <- rep.int(seq_along(object@gene2Nodes), lengths(object@gene2Nodes))[inConditions]
            genes <- split(names(object@gene2Nodes)[geneIndex],
                           factor(conditions[inConditions], levels = unique(whichGO)))
            genes[[.rootCondition(object@parentMapping)]] <- names(object@gene2Nodes)
            lapply(genes[whichGO], sort)
          })

#' @rdname topAnat
#' @export
setMethod("genesInTerm", signature(object = "topAnatPropagatedData", whichGO = "missing"),
          function(object, whichGO) {
            genesInTerm(object, rownames(object@counts))
          })

## The foreground genes with expression data, as sigGenes of topGO
#' @rdname topAnat
#' @export
setMethod("sigGenes", "topAnatPropagatedData", function(object) {
  object@allGenes[object@selected & object@feasible]
})

## The graph of the tested conditions, built as topGO builds it, without the genes of the conditions
#' @rdname topAnat
#' @export
setMethod("graph", "topAnatPropagatedData", function(object) {
  conditions <- sort(unique(unlist(object@gene2Nodes, use.names = FALSE)))
  graph::subGraph(rownames(object@counts), .buildGraph.topology(conditions, object@parentMapping))
})

#' @rdname topAnat
#' @export
setMethod("show", "topAnatPropagatedData", function(object) {
  cat("\n------------------------- topAnatPropagatedData object -------------------------\n")
  cat("\n Description:\n")
  cat("   - ", object@description, "\n")
  cat("\n Ontology:\n")
  cat("   - ", object@ontology, "\n")
  cat("\n", length(object@allGenes), "available genes (all genes of the gene list):\n")
  cat("   - symbol: ", head(object@allGenes, 5), " ...\n")
  cat("   -", sum(object@selected), " significant genes. \n")
  feasibleGenes <- object@allGenes[object@feasible]
  cat("\n", length(feasibleGenes), "feasible genes (genes that can be used in the analysis):\n")
  cat("   - symbol: ", head(feasibleGenes, 5), " ...\n")
  cat("   -", length(sigGenes(object)), " significant genes. \n")
  cat("\n Conditions tested (with at least ", object@nodeSize, " genes):\n")
  cat("   - number of conditions =", nrow(object@counts), "\n")
  cat("\n------------------------- topAnatPropagatedData object -------------------------\n\n")
})
