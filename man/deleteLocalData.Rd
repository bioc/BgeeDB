% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deleteLocalData.R
\name{deleteLocalData}
\alias{deleteLocalData}
\title{Delete local data for the species of the reference class Bgee object}
\usage{
deleteLocalData(myBgeeObject, allDataTypes = FALSE)
}
\arguments{
\item{myBgeeObject}{A Reference Class Bgee object, notably specifying the targeted species}

\item{allDataTypes}{A boolean defining rather data of all datatypes from the selected species
should be deleted. If TRUE, the local sqlite database file is deleted. If FALSE, data
corresponding to the datatype of the \code{myBgeeObject} object will be deleted from the local
sqlite database.}
}
\description{
This function delete data present in the local database. It can delete data
rather from one datatype or all datatypes of one species.
}
\examples{
{
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data <- getData(bgee, experimentId = "SRP007359")
  deleteLocalData(bgee, allDataTypes = TRUE)
}

}
\author{
Julien Wollbrett
}
