% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_fusion_peptide.R
\name{get_fusion_peptide}
\alias{get_fusion_peptide}
\title{translate the fusion sequence to peptides}
\usage{
get_fusion_peptide(obj, n.cores = parallel::detectCores()/2)
}
\arguments{
\item{obj}{a FusionSet object}

\item{n.cores}{number of parallel cores}
}
\value{
a FusionSet object with slot peptide
}
\description{
prepare the short peptides to predict the binding affinity to MHC molecles and determine frame-shift peptide according to 3'-end host gene
}
\examples{
\dontrun{
obj = get_fusion_peptide(obj)
}
}
