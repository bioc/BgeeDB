#' @title Produces an object allowing to perform GO-like enrichment of anatomical terms
#'
#' @description This function produces a topAnatObject, ready to use for gene set enrichment testing. On the expression calls observed in each condition, it is an object of the topGO package, that uses the Uberon ontology instead of the GO ontology. On the expression calls propagated by Bgee, from Bgee 16, BgeeDB runs the classic Fisher test on it itself, without topGO.
#'
#' @details To perform the enrichment test for expression in anatomical structures for each term of Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}), the data are formatted to use the topGO package for testing. This package is interesting because it propagates the mapping of gene to terms to parent terms, and it possesses a pannel of enrichment tests and decorrelation methods. Expert users should be able to use information from the topAnatObject to test enrichment with other packages than topGO.
#' From Bgee 16, the expression calls are either propagated by Bgee or observed in each condition (see the \code{callFileType} argument of \code{loadTopAnatData}). The propagated calls of a condition already hold the genes expressed in its substructures and cell types: \code{runTest} computes on them the classic Fisher test, without decorrelation, for all the conditions at once, without topGO, with the same counts and p-values as the classic algorithm of topGO; the other algorithms and test statistics are refused. \code{termStat}, \code{genesInTerm}, \code{sigGenes} and \code{graph} work on this object as on an object of topGO, as well as \code{makeTable} and \code{printTopAnatGraph}. topGO propagates the observed calls, for its decorrelation algorithms: the classic algorithm is refused on them. A message tells which applies.
#'
#' @param object A topAnatObject, produced by the function \code{topAnat}.
#'
#' @param algorithm The algorithm of topGO to run, see the function \code{runTest} of topGO. Only "classic" on the expression calls propagated by Bgee.
#'
#' @param statistic The test statistic of topGO to use, see the function \code{runTest} of topGO. Only "fisher" on the expression calls propagated by Bgee.
#'
#' @param whichGO The IDs of the conditions, all the tested conditions when missing.
#'
#' @param topAnatData a list including a gene2anatomy list, an organ.relationships list and an organ.names data.frame, produced by the function loadTopAnatData().
#'
#' @param geneList Vector indicating foreground and background genes. Names of the vector indicate the background genes. Values are 1 (gene in foreground) or 0 (gene not in foreground).
#'
#' @param nodeSize Minimum number of genes mapped to a node for it to be tested. Default is 10.
#'
#' @param ... Additional parameters as passed to build topGOdata object in topGO package, for the observed calls.
#'
#' @return topAnatObject: on the expression calls propagated by Bgee, a topAnatPropagatedData class object, ready for the classic Fisher test of \code{runTest}; on the observed calls, a topAnatData class object, ready for gene set enrichment testing with topGO.
#'
#' @author Julien Roux
#'
#' @examples{
#' bgee <- Bgee$new(species="Bos_taurus", dataType="rna_seq")
#' myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000092")
#' geneList <- as.factor(c(rep(0, times=85), rep(1, times=15)))
#' names(geneList) <- c("ENSBTAG00000000011","ENSBTAG00000000014","ENSBTAG00000000016",
#'                      "ENSBTAG00000000026","ENSBTAG00000000039","ENSBTAG00000000040",
#'                      "ENSBTAG00000000042","ENSBTAG00000000050","ENSBTAG00000000056",
#'                      "ENSBTAG00000000064","ENSBTAG00000000067","ENSBTAG00000000071",
#'                      "ENSBTAG00000000072","ENSBTAG00000000080","ENSBTAG00000000081",
#'                      "ENSBTAG00000000084","ENSBTAG00000000091","ENSBTAG00000000099",
#'                      "ENSBTAG00000000111","ENSBTAG00000000123","ENSBTAG00000000132",
#'                      "ENSBTAG00000000153","ENSBTAG00000000162","ENSBTAG00000000163",
#'                      "ENSBTAG00000000169","ENSBTAG00000000179","ENSBTAG00000000197",
#'                      "ENSBTAG00000000199","ENSBTAG00000000202","ENSBTAG00000000203",
#'                      "ENSBTAG00000000204","ENSBTAG00000000213","ENSBTAG00000000215",
#'                      "ENSBTAG00000000223","ENSBTAG00000000224","ENSBTAG00000000225",
#'                      "ENSBTAG00000000236","ENSBTAG00000000250","ENSBTAG00000000251",
#'                      "ENSBTAG00000000252","ENSBTAG00000000253","ENSBTAG00000000261",
#'                      "ENSBTAG00000000274","ENSBTAG00000000277","ENSBTAG00000000279",
#'                      "ENSBTAG00000000285","ENSBTAG00000000286","ENSBTAG00000000287",
#'                      "ENSBTAG00000000289","ENSBTAG00000000297","ENSBTAG00000000305",
#'                      "ENSBTAG00000000312","ENSBTAG00000000328","ENSBTAG00000000335",
#'                      "ENSBTAG00000000341","ENSBTAG00000000343","ENSBTAG00000000354",
#'                      "ENSBTAG00000000355","ENSBTAG00000000356","ENSBTAG00000000365",
#'                      "ENSBTAG00000000372","ENSBTAG00000000379","ENSBTAG00000000380",
#'                      "ENSBTAG00000000382","ENSBTAG00000000396","ENSBTAG00000000404",
#'                      "ENSBTAG00000000405","ENSBTAG00000000406","ENSBTAG00000000411",
#'                      "ENSBTAG00000000425","ENSBTAG00000000434","ENSBTAG00000000435",
#'                      "ENSBTAG00000000438","ENSBTAG00000000448","ENSBTAG00000000451",
#'                      "ENSBTAG00000000454","ENSBTAG00000000456","ENSBTAG00000000457",
#'                      "ENSBTAG00000000459","ENSBTAG00000000462","ENSBTAG00000000469",
#'                      "ENSBTAG00000000470","ENSBTAG00000000484","ENSBTAG00000000497",
#'                      "ENSBTAG00000000501","ENSBTAG00000009707","ENSBTAG00000026266",
#'                      "ENSBTAG00000021992","ENSBTAG00000005353","ENSBTAG00000005333",
#'                      "ENSBTAG00000006424","ENSBTAG00000026972","ENSBTAG00000010799",
#'                      "ENSBTAG00000010799","ENSBTAG00000014614","ENSBTAG00000014614",
#'                      "ENSBTAG00000045757","ENSBTAG00000046332","ENSBTAG00000046332",
#'                      "ENSBTAG00000008394")
#' myTopAnatObject <-  topAnat(myTopAnatData, geneList)
#' }
#'
#' @import topGO graph
#' @export

topAnat <- function(topAnatData, geneList, nodeSize = 10, ... ){

  ## Test if topAnatData not empty
  cat("\nChecking topAnatData object.............\n")
  if( length(topAnatData$gene2anatomy) == 0 ) {
    stop("ERROR: the gene2anatomy list of your topAnatData object is empty.")
  }
  if( length(topAnatData$organ.relationships) == 0 ) {
    stop("ERROR: the organ.relationships list of your topAnatData object is empty.")
  }
  if( length(topAnatData$organ.names[,1]) == 0 ) {
    stop("ERROR: the organ.names data frame of your topAnatData object is empty.")
  }
  ## The kind of expression calls: those of the Bgee releases prior to 16 are observed. The data without
  ## a Bgee object come from the files of the topAnat web page, from Bgee 16 or later.
  callFileType <- if ( is.null(topAnatData$call.file.type) ) "observed" else topAnatData$call.file.type
  propagatedCallsAvailable <- is.null(topAnatData$bgee.object$release) ||
    compareVersion(gsub("_", ".", topAnatData$bgee.object$release), "16.0") >= 0
  if ( callFileType == "propagated" ){
    message("NOTE: the expression calls are propagated by Bgee: BgeeDB runs the classic Fisher test on them, ",
            "without decorrelation, itself, without topGO that would propagate them again.")
  } else if ( propagatedCallsAvailable ){
    message("NOTE: the expression calls are observed in each condition, topGO propagates them to the parent ",
            "conditions: only its decorrelation algorithms can be run on them, the classic algorithm is run ",
            "on the propagated calls.")
  } else {
    message("NOTE: the expression calls are observed in each condition, topGO propagates them to the parent conditions.")
  }

  ## Test if gene list is fine
  cat("\nChecking gene list......................\n")
  if (!is.factor(geneList)){
    geneList <- as.factor(geneList)
  }
  if (length(geneList) == 0){
    stop("ERROR: the gene list provided is empty.")
  }
  if (length(levels(geneList)) != 2){
    stop("ERROR: the gene list provided is not in the right format (should be a named vector including only 0 and 1 values).")
  }
  if (length(geneList) < 100) {
    cat("\nWARNING: Given the low number of genes provided, it is very unlikely that the test will have enough power.\n")
  }
  ## If geneList includes genes not present in topAnatData$gene2anatomy, restrict to these genes
  if (sum(names(geneList) %in% names(topAnatData$gene2anatomy)) != length(geneList)){
    cat("\nWARNING: Some genes in your gene list have no expression data in Bgee, and will not be included in the analysis.", sum(names(geneList) %in% names(topAnatData$gene2anatomy)), "genes in background will be kept.\n")
  }

  if (nodeSize == 0){
    stop("ERROR: the node size parameter has to be at least 1.")
  }
  if (!is.numeric(nodeSize)){
    stop("ERROR: the node size parameter should be a numeric, integer value.")
  }
  nodeSize <- as.integer(nodeSize)

  if ( callFileType == "propagated" ){
    ## This also reports to the user how many genes are in the background / foreground
    topAnatObject <- .topAnatPropagatedData(
                         description = "topAnatPropagatedData class object, ready for the classic Fisher test on the expression calls propagated by Bgee",
                         ontology = "UBERON ontology describing animal anatomical structures",
                         geneList = geneList,
                         nodeSize = nodeSize,
                         parentMapping = topAnatData$organ.relationships,
                         gene2Nodes = topAnatData$gene2anatomy
                         )
    ## The analysis is identified by the gene list, the node size and the numbers of genes of the conditions,
    ## quick to hash, rather than by the expression calls of all the genes
    myAnalysisInfo <- list(geneList, nodeSize, topAnatObject@counts)
  } else {
    ## Building the modified topAnatData object. This also reports to the user how many genes are in the background / foreground
    topAnatObject <- new("topAnatData",
                           description = "topAnatData class object, ready for anatomical ontology enrichment test",
                           ontology = "UBERON ontology describing animal anatomical structures",
                           allGenes = geneList,
                           nodeSize = nodeSize,
                           parentMapping = topAnatData$organ.relationships,
                           gene2Nodes = topAnatData$gene2anatomy
                           )
    topAnatObject@propagatedCallsAvailable <- propagatedCallsAvailable

    ## Create a hash from the topAnatObject and the topAnatData objects
    myAnalysisInfo <- c(topAnatObject, topAnatData)
  }
  ## Use library digest to create a SHA512 hash, that will be used as analysis Id
  analysisId <- digest(myAnalysisInfo, algo = "sha512")
  ## cat(paste0("\nAnalysis Id hash built: ", analysisId, "\n"))

  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## The analysis id is used for our internal statistics. It is a secure hash that does not allow the
  ## identification of the gene lists and datasets used by the user, but will inform us on the
  ## number of distinct analyses that are run with the package.
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ## Send a query to our webservice for statistics purposes of topAnat use (even when cached files are used)
  if (isTRUE(topAnatData$bgee.object$sendStats)){
    myUrl <- paste0(topAnatData$bgee.object$topAnatUrl, "?page=stats&action=launch_top_anat_analysis&api_key=", topAnatData$bgee.object$apiKey, "&analysis_id=", analysisId, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    ## Send query, but no need to wait for the answer
    try(getURL(myUrl, followLocation = TRUE, .opts = list(timeout = 1)), silent=TRUE)
  }

  return(topAnatObject)
}
