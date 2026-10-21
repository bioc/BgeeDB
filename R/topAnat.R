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
#' @param geneList Vector indicating foreground and background genes. Names of the vector indicate the background genes. Values are 1 (gene in foreground) or 0 (gene not in foreground). Either geneList or foregroundGenes must be provided.
#'
#' @param nodeSize Minimum number of genes mapped to a node for it to be tested. Default is 10.
#'
#' @param ... Additional parameters as passed to build topGOdata object in topGO package, for the observed calls.
#'
#' @param background The background of the analysis. With geneList, the gene list is restricted to its genes, the foreground genes outside of it being left out; with foregroundGenes, its genes with expression data are compared to the genes of interest. From Bgee 16, "protein_coding" for the protein coding genes, "polyA" for the genes of the biotypes captured by polyA RNA-Seq, or "observed" for the genes with at least one observation of expression, downloaded from Bgee; or a vector of gene IDs, for a custom background. Default is NULL: "protein_coding" from Bgee 16; for the previous releases and for data built without a Bgee object, no restriction of geneList, and all the genes with expression data for foregroundGenes.
#'
#' @param foregroundGenes The IDs of the genes of interest, to provide instead of geneList: the gene list is then built from the background, its genes of interest in the foreground. Either geneList or foregroundGenes must be provided.
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

topAnat <- function(topAnatData, geneList = NULL, nodeSize = 10, ..., background = NULL, foregroundGenes = NULL){

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
  geneList <- .buildTopAnatGeneList(topAnatData, geneList, foregroundGenes, background, propagatedCallsAvailable)
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

#' Build the gene list of a topAnat analysis
#'
#' Builds the gene list of topGO: a factor with the levels "0" and "1", named by the IDs of the genes of the universe
#' of the analysis, 1 marking the foreground genes. It is built in one of two ways:
#' \itemize{
#'   \item from the gene list given by the user (\code{geneList}): its names are the universe, restricted to the genes
#'   of the background when there is one, the foreground genes outside of it being left out;
#'   \item from the genes of interest given by the user (\code{foregroundGenes}): the universe is then the genes of the
#'   background with expression data in \code{topAnatData}, or all the genes with expression data when there is no
#'   background, the genes of interest in the foreground.
#' }
#' The background is resolved by \code{.topAnatBackgroundGenes}. A message tells the background, and how many genes
#' and genes of interest are kept.
#'
#' @param topAnatData The list produced by \code{loadTopAnatData}, or built from local files by
#' \code{loadTopAnatDataFromFiles}.
#' @param geneList NULL, or the gene list given by the user: a vector of 0 and 1 named by the IDs of the genes.
#' @param foregroundGenes NULL, or the IDs of the genes of interest given by the user. Exactly one of \code{geneList}
#' and \code{foregroundGenes} must be given.
#' @param background The \code{background} argument of \code{topAnat}: NULL, a type of background, or gene IDs.
#' @param propagatedCallsAvailable Whether the release of \code{topAnatData} is Bgee 16 or later.
#' @return The gene list: a factor with the levels "0" and "1", named by gene IDs.
#' @noRd
.buildTopAnatGeneList <- function(topAnatData, geneList, foregroundGenes, background, propagatedCallsAvailable) {
  if ( is.null(geneList) == is.null(foregroundGenes) ){
    stop("ERROR: please provide either a gene list (geneList) or the genes of interest (foregroundGenes), not both.")
  }
  backgroundGenes <- .topAnatBackgroundGenes(topAnatData, background, propagatedCallsAvailable)

  ## From the genes of interest: the universe is the background, restricted to the genes with expression data
  if ( !is.null(foregroundGenes) ){
    foregroundGenes <- unique(as.character(foregroundGenes))
    if ( length(foregroundGenes) == 0 ){
      stop("ERROR: the list of genes of interest provided is empty.")
    }
    universe <- names(topAnatData$gene2anatomy)
    if ( !is.null(backgroundGenes) ){
      universe <- intersect(universe, backgroundGenes$genes)
    }
    inUniverse <- foregroundGenes %in% universe
    message("NOTE: the background is ", if ( is.null(backgroundGenes) ) "all the genes" else backgroundGenes$name,
            " with expression data: ", length(universe), " genes, with ", sum(inUniverse), " of the ",
            length(foregroundGenes), " genes of interest.")
    if ( !any(inUniverse) ){
      stop("ERROR: none of the genes of interest is in the background with expression data.")
    }
    geneList <- factor(as.integer(universe %in% foregroundGenes), levels = c(0, 1))
    names(geneList) <- universe
    return(geneList)
  }

  ## From the gene list of the user: its names are the universe, restricted to the background
  if ( !is.factor(geneList) ){
    geneList <- as.factor(geneList)
  }
  if ( length(geneList) == 0 ){
    stop("ERROR: the gene list provided is empty.")
  }
  if ( length(levels(geneList)) != 2 || is.null(names(geneList)) ){
    stop("ERROR: the gene list provided is not in the right format (should be a named vector including only 0 and 1 values).")
  }
  if ( !is.null(backgroundGenes) ){
    inBackground <- names(geneList) %in% backgroundGenes$genes
    foreground <- as.logical(as.integer(levels(geneList)))[geneList]
    message("NOTE: the background is ", backgroundGenes$name, ": ", sum(inBackground), " of the ", length(geneList),
            " genes of the gene list are kept, ", sum(foreground & !inBackground), " of its ", sum(foreground),
            " foreground genes are left out.")
    if ( !any(foreground & inBackground) ){
      stop("ERROR: none of the foreground genes is in the background.")
    }
    geneList <- geneList[inBackground]
  }
  geneList
}

#' Resolve the background of a topAnat analysis
#'
#' Resolves the \code{background} argument of \code{topAnat} into the IDs of the genes of the background:
#' \itemize{
#'   \item NULL: the "protein_coding" background from Bgee 16 with a Bgee object; otherwise no background, for the
#'   previous releases and for data built from local files without a Bgee object;
#'   \item "protein_coding", "polyA" or "observed": downloaded from Bgee by \code{.loadTopAnatBackground}, which needs
#'   a Bgee object and Bgee 16 or later;
#'   \item any other vector: the IDs of the genes of a custom background.
#' }
#'
#' @param topAnatData The list produced by \code{loadTopAnatData}, or built from local files by
#' \code{loadTopAnatDataFromFiles}.
#' @param background The \code{background} argument of \code{topAnat}.
#' @param propagatedCallsAvailable Whether the release of \code{topAnatData} is Bgee 16 or later.
#' @return NULL when there is no background; otherwise a list with \code{genes}, the IDs of the genes of the
#' background, and \code{name}, its description for the messages.
#' @noRd
.topAnatBackgroundGenes <- function(topAnatData, background, propagatedCallsAvailable) {
  if ( is.null(background) && !is.null(topAnatData$bgee.object$release) && propagatedCallsAvailable ){
    background <- "protein_coding"
  }
  if ( is.null(background) ){
    return(NULL)
  }
  if ( length(background) == 1 && background %in% c("protein_coding", "polyA", "observed") ){
    if ( is.null(topAnatData$bgee.object$release) ){
      stop("ERROR: the background \"", background, "\" is downloaded from Bgee, which needs a Bgee object: ",
           "please provide the IDs of the genes of the background.")
    }
    if ( !propagatedCallsAvailable ){
      stop("ERROR: the background \"", background, "\" is available from Bgee 16: ",
           "please provide the IDs of the genes of the background.")
    }
    return(list(genes = .loadTopAnatBackground(topAnatData$bgee.object, background),
                name = paste0("the \"", background, "\" genes of Bgee")))
  }
  list(genes = unique(as.character(background)), name = "the genes provided")
}

## The genes of a background of topAnat, from Bgee 16: downloaded from the Bgee API, and kept in the download
## directory as the other files of topAnat
.loadTopAnatBackground <- function(myBgeeObject, background) {
  backgroundFileName <- paste0("topAnat_Background_", myBgeeObject$speciesId, "_", toupper(background), ".tsv")
  backgroundFile <- file.path(myBgeeObject$pathToData, backgroundFileName)
  if (file.exists(backgroundFile)) {
    cat(paste0("\nNOTE: a background file was found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))
  } else {
    myUrl <- paste0(myBgeeObject$topAnatUrl, "?page=r_package&action=get_background_genes&display_type=tsv&species_list=",
                    myBgeeObject$speciesId, "&background_type=", background, "&api_key=", myBgeeObject$apiKey,
                    "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    cat(paste0("\nRetrieving the genes of the background from Bgee (", myUrl, ")\n"))
    ## A query refused by the server makes download.file warn and fail without saying why: the answer of the server
    ## is then retrieved below
    success <- tryCatch(suppressWarnings(bgee_download_file(url = myUrl, destfile = paste0(backgroundFile, ".tmp"))),
                        error = function(e) 1)
    if (success != 0) {
      unlink(paste0(backgroundFile, ".tmp"))
      serverAnswer <- try(getURL(myUrl), silent = TRUE)
      if (inherits(serverAnswer, "try-error")) {
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      }
      stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n",
                  trimws(serverAnswer)))
    }
    ## The 5 blank lines ending a response tell a complete transfer
    if (!.isCompleteResponse(paste0(backgroundFile, ".tmp"))) {
      file.remove(paste0(backgroundFile, ".tmp"))
      stop(paste0("File ", backgroundFileName, " is truncated, there may be a temporary problem with the Bgee webservice."))
    }
    file.rename(paste0(backgroundFile, ".tmp"), backgroundFile)
  }
  read.table(backgroundFile, header = TRUE, sep = "\t", quote = "", comment.char = "", colClasses = "character")$GENE_ID
}
