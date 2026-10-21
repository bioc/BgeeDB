## The S4 classes of topAnat. The files of the package are collated in alphabetical order: the classes are
## defined in this file, first, before their methods.

######################## topAnatData class ########################
## This class is an extension of the topGOdata class from
## the topGO package. The additional fields are:
## - parentMapping: a list describing parent-child
##   relationships in the ontology to be used in place of
##   the Gene Ontology
## - gene2Nodes: a list providing the mapping of genes to
##   term from the new ontology

## Declare topAnatData as an extension of topGOdata class
setClass("topAnatData",
         representation = representation(
           ## Added slot to give us the structure of the ontology
           parentMapping = "character",
           ## Added slot to give us the mapping of genes to terms of the ontology
           gene2Nodes = "character",
           ## Added slot telling whether the release of the calls has propagated calls, from Bgee 16:
           ## the classic algorithm is then only run on them
           propagatedCallsAvailable = "logical"),
         contains = "topGOdata"
)

######################## topAnatPropagatedData class ########################
## The expression calls propagated by Bgee already hold, in each condition, the genes expressed in its
## substructures and cell types. The classic Fisher test is then computed from the numbers of genes and of
## foreground genes of each condition, for all the conditions at once, without the object of topGO, whose
## environment of the genes of each condition is long to build and large for the species with many conditions.
## The genes are counted as topGO counts them, including those given several times in the gene list.
setClass("topAnatPropagatedData",
         representation = representation(
           description = "character",
           ontology = "character",
           ## The genes of the gene list, whether each is in the foreground, and whether it has expression data
           allGenes = "character",
           selected = "logical",
           feasible = "logical",
           ## Minimum number of genes of a tested condition
           nodeSize = "integer",
           ## The conditions where each gene with expression data is expressed
           gene2Nodes = "list",
           ## The parents of each condition
           parentMapping = "list",
           ## The numbers of genes (Annotated) and of foreground genes (Significant) of each tested condition,
           ## in the order of the IDs of the conditions, which are the row names
           counts = "data.frame")
)
