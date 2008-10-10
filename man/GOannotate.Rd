\name{GOannotate}
\alias{GOannotate}

\title{ Retrieve the GO annotations for a list of gene names }
\description{
  Retrieves the GO annotations for a list of gene names with the help of the Bioconductor package biomaRt. Queries only for Homo sapiens.
}
\usage{
GOannotate(vecOfGeneNames, whichOnto)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vecOfGeneNames}{ a vector of gene names }
  \item{whichOnto}{ One of the three GO ontologies: \code{ "biological_process"}, \code{ "molecular_function"} or \code{ "cellular_component"} }
}

\value{
  A named list of GO terms of the given ontology.
}
\references{ http://geneontology.org/}

\seealso{ \code{\link{gseaAnalysis}} }

\keyword{ internal }
