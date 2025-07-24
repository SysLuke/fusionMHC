#' This function implements the prediction of binding affinity to MHC molecules for each peptide
#' @description predict the binding of each peptide in obj@peptide to MHC class I molecules, i.e, HLA-A0201 and HLA-A1101 using netMHC-4.0 (https://services.healthtech.dtu.dk/services/NetMHC-4.0/).
#' @param obj a FusionSet object
#' @param path.to.netMHC path to netMHC-4.0
#' @param temp.dir specify temporary directory
#' @param n.cores number of parallel cores
#' @return returns a \code{FusionSet} object with slot occurrence and renew slot peptide with extra columns prefix with HLA-*).
#' the binding details were added to obj@peptide, prefix with corresponding HLA type. occurrence matrix represents the peptides with MHC binding affinity for each HLA type in each sample.
#'
#' @export
#' @examples
#' \dontrun{
#' obj = predict_MHC_affinity(obj, path.to.netMHC = "path/to/netMHC-4.0/")
#' }

predict_MHC_affinity <- function(obj, path.to.netMHC = "path/to/netMHC-4.0/", temp.dir = tempdir(), n.cores = parallel::detectCores()/2){
  pt = obj@peptide
  if(!"Fusion.peptide.MHC" %in% names(pt) | nrow(pt) == 0){
    stop("no valiable records in obj@peptide. run get_fusion_peptide first.")
  }

  id_pep = which(names(pt) == 'Fusion.peptide.MHC')
  id_shift = which(names(pt) == 'Frame.shift')
  out = parallel::mclapply(1:nrow(pt), function(row){
    x = pt[row,]
    pep = x[[id_pep]]

    input.fa = paste0(temp.dir,"/file_",row,".fsa")
    xls_file = paste0(temp.dir,"/file_",row,".xls")
    log_file = paste0(temp.dir,"/file_",row,".log")

    sink(input.fa)
    cat(">",row,"\n",sep = "")
    cat(pep)
    sink()

    cmd = paste0(path.to.netMHC, "/Linux_x86_64/bin/netMHC -a HLA-A0201,HLA-A1101 -l 8,9,10,11 -rth 0.5 -rlt 2 -tdir ",temp.dir, " -version ",path.to.netMHC, "/data/version -hlalist ", path.to.netMHC, "/data/allelelist -thrfmt ",
                 path.to.netMHC, "/data/threshold/%s.thr -syn ", path.to.netMHC, "/data/synlists/%s.synlist -inptype 0 -rdir ", path.to.netMHC, "/Linux_x86_64/ -f ",
                 input.fa," -xls -xlsfile ", xls_file," 1> ", log_file," 2>&1")
    system(cmd)

    df = data.table::fread(xls_file, data.table = F, check.names = T)
    df$n_aa = nchar(df$Peptide)
    df$shift = x[[id_shift]]
    df = df[which((df$shift %in% c("No",'Short.match') & df$Pos+df$n_aa >= 11 & df$Pos<=9) |
                    (df$shift == "Yes" & df$Pos+df$n_aa >= 11)),]
    o <- NULL
    if(nrow(df) > 0){
      if(df$ID[1] == row){
        o0 = data.frame(row = row)

        ix = which.min(df$Rank)
        o1 = data.frame(n_sb = sum(df$Rank <= 0.5), n_wb = sum(df$Rank <= 2 & df$Rank > 0.5),
                        Affinity=df$nM[ix], Rank=df$Rank[ix], aaseq = df$Peptide[ix],
                        n_aa = df$n_aa[ix])
        o1$binding = ifelse(o1$n_sb + o1$n_wb, "Yes", 'No')
        names(o1) = paste0('HLA-A0201:',names(o1))

        ix = which.min(df$Rank.1)
        o2 = data.frame(n_sb = sum(df$Rank.1 <= 0.5), n_wb = sum(df$Rank.1 <= 2 & df$Rank.1 > 0.5),
                        Affinity=df$nM.1[ix],Rank=df$Rank.1[ix], aaseq = df$Peptide[ix],
                        n_aa = df$n_aa[ix])
        o2$binding = ifelse(o2$n_sb + o2$n_wb, "Yes", 'No')
        names(o2) = paste0('HLA-A1101:',names(o2))

        o = cbind(o0,o1,o2)
      }
    }
    system(paste0("rm -f ",input.fa, " ", xls_file, " ", log_file))
    return(o)
  }, mc.cores = n.cores)

  out = do.call(rbind, out)
  pt[out$row,names(out)[-1]] = out[,-1]

  # count the occurrence
  counts = obj@count
  oList = list()
  for(type in c("HLA-A0201","HLA-A1101")){
    binding_peps  = unique(pt$chimericRNA[which(pt[,paste0(type,":binding")]=='Yes')])
    idx = match(binding_peps, rownames(counts))
    data = counts[idx, ]
    data[data > 0] = 1

    oList[[type]] = data
  }

  obj@peptide = pt
  obj@occurrence = oList
  return(obj)
}
