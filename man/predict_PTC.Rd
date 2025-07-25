% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_PTC.R
\name{predict_PTC}
\alias{predict_PTC}
\title{detects the peptides harboring premature termination codons}
\usage{
predict_PTC(obj, gtf.ele.list, fa.set)
}
\arguments{
\item{obj}{a FusionSet object}

\item{gtf.ele.list}{output of parse_gtf}

\item{fa.set}{output of parse_fasta}
}
\value{
a FusionSet object with slot PTC
}
\description{
this function detects the fusion peptides which harbor premature termination codons (PTC) defined by stop codon terminated at last 55 bps of the second to last exon and last exon.
}
\examples{
\dontrun{
gtf.ele.list = parse_gtf(gtf.file='path/to/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz')
fa.set = parse_fasta(fasta.file='path/to/GRCh38.fa')
obj = predict_PTC(obj, gtf.ele.list, fa.set)
}
}
