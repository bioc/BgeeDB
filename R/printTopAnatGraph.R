#' @title Draw the graph of the conditions significant in a topAnat analysis
#'
#' @description This function draws, in a PDF file, the conditions most significant in a topAnat analysis and all their ancestors, with their names, their p-values, and their numbers of significant and annotated genes. The most significant conditions are drawn as boxes, colored by p-value, their ancestors as ellipses.
#'
#' @details The graph is drawn with the topGO package, and needs the Rgraphviz package, to install with \code{BiocManager::install("Rgraphviz")}. It works for the expression calls propagated by Bgee and for the observed ones.
#'
#' @param topAnatData A list produced by the function \code{loadTopAnatData}.
#'
#' @param topAnatObject A topAnatData class object, produced by the function \code{topAnat}.
#'
#' @param results A topGOresult object, produced by the function \code{runTest} of topGO.
#'
#' @param firstSigNodes Number of most significant conditions to draw, with their ancestors. Default is 10.
#'
#' @param fileName Path to the PDF file to write. Default is "topAnat_graph.pdf".
#'
#' @param useInfo The information written in each node: "all" (default), the name, the p-value and the numbers of significant and annotated genes; "np", the name and the p-value; "def", the name; "pval", the p-value; "counts", the numbers of genes; "none", only the ID.
#'
#' @return The path to the PDF file, invisibly.
#'
#' @author Julien Wollbrett
#'
#' @examples{
#' \dontrun{
#'   bgee <- Bgee$new(species = "Danio_rerio", release = "15.2")
#'   myTopAnatData <- loadTopAnatData(bgee)
#'   data(geneList)
#'   myTopAnatObject <- topAnat(myTopAnatData, geneList)
#'   results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')
#'   printTopAnatGraph(myTopAnatData, myTopAnatObject, results, firstSigNodes = 5)
#' }
#' }
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par
#' @export
printTopAnatGraph <- function(topAnatData, topAnatObject, results, firstSigNodes = 10,
                              fileName = "topAnat_graph.pdf",
                              useInfo = c("all", "np", "def", "pval", "counts", "none")){
  if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
    stop("ERROR: the Rgraphviz package is needed to draw the graph, please install it with BiocManager::install(\"Rgraphviz\").")
  }
  useInfo <- match.arg(useInfo)
  pValues <- score(results)
  if (length(pValues) == 0) {
    stop("ERROR: the results object is empty.")
  }

  ## The most significant conditions and all their ancestors, the root at the top
  sigTerms <- sort(pValues)[seq_len(min(firstSigNodes, length(pValues)))]
  dag <- reverseArch(inducedGraph(graph(topAnatObject), names(sigTerms)))
  conditions <- graph::nodes(dag)

  ## The information written in each node, with the names of the conditions where topGO would use
  ## the definitions of the GO terms
  conditionNames <- topAnatData$organ.names[match(conditions, topAnatData$organ.names[, 1]), 2]
  pval <- format.pval(pValues[conditions], digits = 3, eps = 1e-20)
  termCounts <- termStat(topAnatObject, conditions)
  counts <- apply(termCounts[, c("Significant", "Annotated")], 1, paste, collapse = " / ")
  nodeInfo <- switch(useInfo,
                     none = NULL,
                     pval = pval,
                     def = conditionNames,
                     counts = counts,
                     np = paste(conditionNames, pval, sep = "\\\n"),
                     all = paste(conditionNames, pval, counts, sep = "\\\n"))

  ## topGO draws the graph with functions of Rgraphviz it does not import: they must be attached
  if (!"package:Rgraphviz" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("Rgraphviz"))
    on.exit(detach("package:Rgraphviz"), add = TRUE)
  }
  completeDag <- GOplot(dag, sigNodes = pValues[conditions], genNodes = names(sigTerms),
                        wantedNodes = names(sigTerms), showEdges = TRUE, useFullNames = TRUE,
                        oldSigNodes = NULL, nodeInfo = nodeInfo)
  pdf(file = fileName, width = 25, height = 25)
  par(mai = rep(0, 4))
  plot(completeDag)
  dev.off()
  cat(paste0("\nGraph of ", length(conditions), " conditions written in ", fileName, "\n"))
  return(invisible(fileName))
}
