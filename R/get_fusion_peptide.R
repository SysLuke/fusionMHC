#' translate the fusion sequence to peptides
#' @description prepare the short peptides to predict the binding affinity to MHC molecles and determine frame-shift peptide according to 3'-end host gene
#' @param obj a FusionSet object
#' @param n.cores number of parallel cores
#' @return a FusionSet object with slot peptide
#'
#' @export
#' @examples
#' \dontrun{
#' obj = get_fusion_peptide(obj)
#' }

get_fusion_peptide <- function(obj, n.cores = parallel::detectCores()/2){
  win_size = 10

  fd = tibble::rownames_to_column(obj@feature.data,var='chimericRNA')
  id_seq = which(names(fd) == 'Fusion.sequence')
  id_entry5 = which(names(fd) == 'Protein.entry.5')
  id_entry3 = which(names(fd) == 'Protein.entry.3')
  id_type = which(names(fd) == 'Fusion.type')
  id_chimericRNA = which(names(fd) == 'chimericRNA')

  message("conduct frame-shift translation and determine correct frame-shift peptides to predict the binding affinity")
  out = parallel::mclapply(1:nrow(fd), function(row){
    x = fd[row,]
    x1 = gsub("\\*","",x[[id_seq]])
    x2 = gsub("\\*.*","",x[[id_seq]])

    o = data.frame()
    for(i in 0:2){
      m = Biostrings::DNAString(x1, start = i+1)
      pep0 = as.character(suppressWarnings(Biostrings::translate(m, if.fuzzy.codon = 'X')))
      codon_loc = ifelse(stringr::str_detect(pep0,'\\*'),'Yes','No')
      pep = gsub("\\*.*","",pep0)
      fp = floor((nchar(x2)-i)/3)

      max_score = 0
      shift = 'yes'
      if(x[[id_type]] == 'EE'){
        # best match
        ref_i = as.character(ref_pep[x[[id_entry5]]])
        if(Biostrings::width(ref_i) >= win_size & nchar(pep) >= win_size){
          pep_sub = substr(pep,1,win_size)
          for(j in 1:(nchar(ref_i)-win_size+1)){
            ref_sub = substr(ref_i,j,j+9)
            max_score = max(max_score,
                            sum(stringr::str_split(pep_sub,"")[[1]] == stringr::str_split(ref_sub,"")[[1]]), na.rm = T)
            if(max_score >= 10){
              break
            }
          }
        }
      }
      o = rbind(o, data.frame(chimericRNA = x[[id_chimericRNA]], Frame.num = i, Full.peptide = pep0, Effective.peptide = pep, Stop.codon = codon_loc, Stop.pos = nchar(pep),
                              Fusion.point = fp, Match.score = max_score))
    }

    o$Cross.fusion = ifelse(o$Stop.pos >= o$Fusion.point, 'Yes','No')
    o$Is.best = ifelse(x[[id_type]] %in% c('non-coding','non-EE'),'ND','No')
    o[which.max(o$Match.score),"Is.best"] = ifelse(x[[id_type]] == 'EE' & max(o$Match.score) >= 6,'Yes','No')
    o$Frame.shift = 'Yes'
    o[o$Stop.pos - o$Fusion.point <= 7,'Frame.shift'] = 'Short.match'
    o = o[which(o$Cross.fusion == 'Yes' & o$Is.best %in% c('Yes','ND')),]
    o$Cross.fusion = NULL

    if(nrow(o) > 0){
      for(i in 1:nrow(o)){
        if(!is.na(x[[id_entry3]]) & o$Stop.pos[i] - o$Fusion.point[i] > 7){
          ref_i = as.character(ref_pep[x[[id_entry3]]])
          seq_sub3 = substr(o$Effective.peptide[i],o$Fusion.point[i]+1,nchar(o$Effective.peptide[i]))
          o[i,'Frame.shift'] = ifelse(stringr::str_detect(ref_i, seq_sub3), 'No', 'Yes')
        }
      }
      o$Fusion.peptide.MHC = ifelse(o$Frame.shift == 'No',
                                    substr(o$Effective.peptide, o$Fusion.point-9,o$Fusion.point+1+9),
                                    substr(o$Effective.peptide, o$Fusion.point-9,nchar(o$Effective.peptide)))
    }
    return(o)
  }, mc.cores = n.cores)

  out = do.call(rbind, out)
  rownames(out) = NULL

  obj@peptide = out
  return(obj)
}
