#' prepare the reference sequence
#' @description prepare the reference sequence
#' @param fasta.file fasta file with corresponding version referring to FusionCatcher output
#' @return a DNAStringSet object
#'
#' @export
#' @examples
#' \dontrun{
#' fa.set = parse_fasta(fasta.file='path/to/GRCh38.fa')
#' }

parse_fasta <- function(fasta.file){
  x <- Biostrings::readDNAStringSet(fasta.file, format = "fasta")
  x@ranges@NAMES = gsub(" .*","",x@ranges@NAMES)
  return(x)
}
