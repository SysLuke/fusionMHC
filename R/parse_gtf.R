#' prepare the gene annotation from the gtf file
#' @description extract the exon elements from the gtf file
#' @param gtf.file gtf file referring to version used in FusionCatcher
#' @return a gene annotation list including exon elements
#'
#' @export
#' @examples
#' \dontrun{
#' # downloaded from //ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz.
#' # you should download the corresponding version referring to FusionCatcher output
#' gtf.ele.list = parse_gtf(gtf.file='path/to/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz')
#' }

parse_gtf <- function(gtf.file){
  gtf = rtracklayer::import(gtf.file)
  gtf_exon = gtf[gtf@elementMetadata$type =='exon' & gtf@elementMetadata$tag == 'MANE Select',]
  gtf_exon$exon_number = as.numeric(gtf_exon$exon_number)
  gtf_exon = as.data.frame(gtf_exon)
  gtf_exon = gtf_exon[,c('start','end','strand','gene_id','exon_number')]
  gtf_exon_full_list = split(gtf_exon,f=gtf_exon$gene_id)
  gtf_exon_full_list = lapply(gtf_exon_full_list,
                              FUN = function(x){
                                x$gene_id = NULL
                                rownames(x) = NULL
                                return(x)})
  gtf_exon_ptc_list = lapply(gtf_exon_full_list, function(x){
    m = max(x$exon_number)
    x = x[x$exon_number >= m-1,]
    if(x$strand[1] == '-'){
      x[x$exon_number == m-1, 'end'] = x[x$exon_number == m-1, 'start'] + 54
    } else{
      x[x$exon_number == m-1, 'start'] = x[x$exon_number == m-1, 'end'] - 54
    }
    x = x[,c('start','end')]
    rownames(x) = NULL
    colnames(x) = NULL
    return(x)
  })
  x = list(full = gtf_exon_full_list, ptc = gtf_exon_ptc_list)
  return(x)
}
