#' detects the peptides harboring premature termination codons
#' @description this function detects the fusion peptides which harbor premature termination codons (PTC) defined by stop codon terminated at last 55 bps of the second to last exon and last exon.
#' @param obj a FusionSet object
#' @param gtf.ele.list output of parse_gtf
#' @param fa.set output of parse_fasta
#' @return a FusionSet object with slot PTC
#'
#' @export
#' @examples
#' \dontrun{
#' gtf.ele.list = parse_gtf(gtf.file='path/to/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz')
#' fa.set = parse_fasta(fasta.file='path/to/GRCh38.fa')
#' obj = predict_PTC(obj, gtf.ele.list, fa.set)
#' }

predict_PTC <- function(obj, gtf.ele.list, fa.set){
  gtf_exon_ptc_list = gtf.ele.list$ptc
  gtf_exon_full_list = gtf.ele.list$full

  fd = obj@feature.data
  fd = fd[which(fd$Predicted.effect %in% c('out-of-frame', 'CDS(truncated)/UTR')),]
  fd = tibble::rownames_to_column(fd, var='chimericRNA')

  pt = obj@peptide
  pt = pt[which(pt$chimericRNA %in% fd$chimericRNA),]

  message("only fusion types: out-of-frame and CDS(truncated)/UTR can be supported.")

  df = merge(fd,pt,by='chimericRNA')
  df$Position.3 = as.numeric(df$Position.3)
  df$PTC = ifelse(df$Hostgene3 %in% names(gtf_exon_ptc_list) & df$Chr.3 %in% names(fa.set), "No", "reference gene or chromosome not found")
  print(table(df$PTC))

  df = df[which(df$PTC == "No"),]
  if(nrow(df) == 0) stop('no records left to predict.')
  message(nrow(df), " records to predict PTC.")

  out = parallel::mclapply(1:nrow(df), function(row){
    x = df[row,]
    if(x$Stop.codon == 'Yes'){ # there's stop codon found in 3-end fusion sequence, go to predict PTC
      base_seq = substr(x$Fusion.sequence,x$Frame.num+1,3*x$Stop.pos+1)
      len = nchar(gsub(".*\\*","",base_seq))
      pos2 = ifelse(x$Strand.3 == "+", x$Position.3+len+2, x$Position.3-len-2)
      z = gtf_exon_ptc_list[[x$Hostgene3]]
      n = sum(z[,1] <= pos2 & z[,2] >= pos2)
      if(n>0) x$PTC='Yes'
    } else{ # no stop codon found in 3-end fusion sequence, continue to translate until first stop codon, then predict PTC
      ref_seq = fa.set[[x$Chr.3]]
      y = gtf_exon_full_list[[x$Hostgene3]]
      y = y[order(y$exon_number),]
      idx = y$exon_number[which(y$start<=x$Position.3 & y$end>=x$Position.3)]
      y = y[y$exon_number >= idx,]
      if(nrow(y) == 0){
        x$PTC = "outside the gene"
      } else{
        pre_seq = substr(gsub("\\*.*",'\\*',x$Fusion.sequence),x$Frame.num+1,1000)
        if(x$Strand.3 == "+"){
          y$start[1] = x$Position.3
        } else{
          y$end[1] = x$Position.3
        }
        for(i in 1:nrow(y)){
          suf_seq = Biostrings::subseq(ref_seq,start = y$start[i],end=y$end[i])
          if(x$Strand.3 == "-"){
            suf_seq = Biostrings::reverseComplement(suf_seq)
          }
          pre_seq = paste0(pre_seq,suf_seq)
          pre_seq2 = Biostrings::DNAString(x=gsub("\\*","",pre_seq))
          tmp_aa = as.character(suppressWarnings(Biostrings::translate(pre_seq2, if.fuzzy.codon = 'X')))
          if(stringr::str_detect(tmp_aa, "\\*")){
            x$Fusion.sequence = as.character(pre_seq)
            x$Full.peptide = tmp_aa
            x$Effective.peptide = gsub("\\*.*","",tmp_aa)
            x$Stop.codon = 'Yes'
            x$Stop.pos = nchar(x$Effective.peptide)
            x$Fusion.seqLen = nchar(pre_seq)-1

            base_seq = substr(x$Fusion.sequence,x$Frame.num+1,3*x$Stop.pos+1)
            len = nchar(gsub(".*\\*","",base_seq))
            pos2 = ifelse(x$Strand.3 == "+", x$Position.3+len+2, x$Position.3-len-2)

            z = gtf_exon_ptc_list[[x$Hostgene3]]
            n = sum(z[,1] <= pos2 & z[,2] >= pos2)
            if(n>0) x$PTC='Yes'
            break
          }
        }
      }
    }
    return(x)
  })
  out = do.call(rbind, out)

  out = out[which(out$PTC == "Yes"),]
  out$PCT = NULL
  message(nrow(out), " PTC events found in total.")
  obj@PTC = out
  return(obj)
}
