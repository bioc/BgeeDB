## The methods of the topAnatData class, defined in AllClasses.R: the object of topGO on the observed expression
## calls, that topGO propagates. The graph of the conditions built here also serves the topAnatPropagatedData class.

#################################################################
## Sets of functions written with help of Adrian Alexa (author of the topGO package, personnal communication) to make topGO work with another ontology than the Gene Ontology

.annFUN.gene2Nodes <- function(feasibleGenes = NULL, gene2Nodes) {
  ## Restrict the mappings to the feasibleGenes set
  if(!is.null(feasibleGenes))
    gene2Nodes <- gene2Nodes[intersect(names(gene2Nodes), feasibleGenes)]

  ## Throw-up the genes which have no annotation
  if(any(is.na(gene2Nodes)))
    gene2Nodes <- gene2Nodes[!is.na(gene2Nodes)]

  gene2Nodes <- gene2Nodes[sapply(gene2Nodes, length) > 0]

  ## Get all the Terms and keep a one-to-one mapping with the genes
  allNodes <- unlist(gene2Nodes, use.names = FALSE)
  geneID <- rep(names(gene2Nodes), sapply(gene2Nodes, length))

  return(split(geneID, allNodes))
}
#################################################################

.buildGraph.topology <- function(knownNodes, parentMapping) {
  ## first build the lookUp table for the terms
  nodeLookUp <- new.env(hash = TRUE, parent = emptyenv())

  ## warping functions for a easier acces to the lookUp table
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLookUp)
  }

  ## get the root of the ontology
  GENE.ONTO.ROOT <- setdiff(unique(unlist(parentMapping, use.names = FALSE)),
      names(parentMapping))
  if (length(GENE.ONTO.ROOT) > 1){
    cat('\n Multiple roots were found... Possibly some nodes are not listed in the organ.names data frame\n')
  }

  ## we read all the database once and the access the list
  adjLookUP <- as.list(parentMapping)

  ## we use an environment of environments to store edges: (this way is faster)
  ## in the end we will coerce it to a list of list and build a graphNEL obj.
  edgeEnv <- new.env(hash = TRUE, parent = emptyenv())

  ## add the arc (u --> v) to edgeEnv of type :
  envAddEdge <- function(u, v) {
    assign(v, 0, envir = get(u, envir = edgeEnv))
  }

  ## recursivly build the induced graph starting from one node
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(isNodeInDAG(node))
      return(1)

    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = TRUE, parent = emptyenv()), envir = edgeEnv) # adj list

    if(node == GENE.ONTO.ROOT)
      return(2)

    adjNodes <- adjLookUP[[node]]

    ## debuging point! should not happen!
    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')

    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x)
      buildInducedGraph(x)
    }

    return(0)
  }

  ## we start from the most specific nodes
  lapply(knownNodes, buildInducedGraph)

  ## now we must transform our env into a Graph structure
  ## for now we use lapply, later we can do it with eapply
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })

  ## now we can build the graphNEL object
  return(new('graphNEL',
             nodes = .graphNodes,
             edgeL = .edgeList,
             edgemode = 'directed'))
}
######################## topAnatData class, defined in AllClasses.R ########################

setMethod("initialize", "topAnatData",
          function(.Object,
                   ## which Ontology to be used
                   ontology = "Uberon",
                   ## a named numeric or factor, the names are the genes ID
                   allGenes,
                   ## function to select the signif. genes
                   geneSelectionFun = NULL,
                   ## description of the class
                   description = character(0),
                   ## minimum node size
                   nodeSize = 1,
                   ## the child-parent relationship
                   parentMapping,
                   ## annotation data
                   gene2Nodes,
                   ## additional parameters
                   ...) {

            .Object@description <- description
            .Object@ontology <- ontology

            ## some checking
            if(is.null(names(allGenes)))
              stop("allGenes must be a named vector")

            if(!is.factor(allGenes) && !is.numeric(allGenes))
              stop("allGenes should be a factor or a numeric vector")

            .Object@allGenes <- names(allGenes)

            if(is.factor(allGenes)) {
              if(length(levels(allGenes)) != 2)
                stop("allGenes must be a factor with 2 levels")
              .Object@allScores <- factor(as.character(allGenes))
              .Object@geneSelectionFun <- function(x) {
                return(as.logical(as.integer(levels(x)))[x])
              }
            }
            else {
              .Object@allScores <- as.numeric(allGenes)

              ## function to select which genes are significant
              if(is.null(geneSelectionFun))
                warning("No function to select the significant genes provided!")
              .Object@geneSelectionFun <- geneSelectionFun
            }

            ## size of the nodes which will be pruned
            .Object@nodeSize = as.integer(max(nodeSize, 1))

            ## this function is returning a list of terms from the specified ontology
            ## with each entry being a vector of genes
            cat("\nBuilding most specific Ontology terms... ")
            mostSpecificTerms <- .annFUN.gene2Nodes(feasibleGenes = .Object@allGenes, gene2Nodes = gene2Nodes)
            cat(" ( ", length(mostSpecificTerms), " Ontology terms found. )\n")

            ## the graph is build starting from the most specific terms
            cat("\nBuilding DAG topology................... ")
            g <- .buildGraph.topology(names(mostSpecificTerms), parentMapping)
            cat(" ( ",  graph::numNodes(g), " Ontology terms and ", graph::numEdges(g), " relations. )\n")

            ## probably is good to store the leves but for the moment we don't
            .nodeLevel <- buildLevels(g, leafs2root = TRUE)

            ## annotate the nodes in the graph with genes
            cat("\nAnnotating nodes (Can be long).......... ")
            g <- topGO:::mapGenes2GOgraph(g, mostSpecificTerms, nodeLevel = .nodeLevel) ## leafs2root

            ## select the feasible genes
            gRoot <- getGraphRoot(g)
            feasibleGenes <- ls(graph::nodeData(g, n = gRoot, attr = "genes")[[gRoot]])
            cat(" ( ", length(feasibleGenes), " genes annotated to the Ontology terms. )\n")

            .Object@feasible <- .Object@allGenes %in% feasibleGenes


            ## prune the GO graph
            if(.Object@nodeSize > 1) {
              cc <- .countsInNode(g, graph::nodes(g))
              .Object@graph <- graph::subGraph(names(cc)[cc >= .Object@nodeSize], g)
            } else {
              .Object@graph <-  g
            }

            .Object
          })

## Runs a test of topGO on the observed expression calls, that topGO propagates, after checking that it suits
## them: from Bgee 16, only the decorrelation algorithms, the classic algorithm being run on the calls
## propagated by Bgee.
#' @rdname topAnat
#' @export
setMethod("runTest", signature(object = "topAnatData", algorithm = "character", statistic = "character"),
          function(object, algorithm, statistic, ...) {
            if (algorithm == "classic" && isTRUE(object@propagatedCallsAvailable)) {
              stop(paste0("ERROR: from Bgee 16, the classic algorithm, without decorrelation, is run on the expression ",
                          "calls propagated by Bgee: please retrieve them with loadTopAnatData(..., callFileType = ",
                          "\"propagated\"). The observed calls are for the decorrelation algorithms."))
            }
            message("Running the ", algorithm, " algorithm of topGO on the observed expression calls, ",
                    "propagated by topGO.")
            callNextMethod()
          })
