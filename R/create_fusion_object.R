#' create a fusion object from raw reads count
#' @description import the FusionCatcher output to create a FusionSet object
#' @param input.dir directory to all final-list_candidate-fusion-genes.txt renamed to sample_name.txt (character)
#' @param min.count filter fusion events at least min.count spanning unique reads (numeric)
#' @param info sample information with sample id as row names (data frame)
#' @return a FusionSet object.
#'
#' @import dplyr
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' input.dir = system.file("extdata","fusioncatcher_output", package = "fusionMHC")
#' info = data.frame(row.names = c("sample_1_CA","sample_1_NAT","sample_2_CA",
#'                  "sample_2_NAT","sample_3_CA","sample_3_NAT"),
#'                   group = rep(c('CA','NAT'),3))
#' obj = create_fusion_object(input.dir, info=info)
#' }

create_fusion_object <- function(input.dir, min.count = 3, info = NULL){
  fs = list.files(input.dir, full.names = T, include.dirs = T)
  fList = lapply(fs, function(f){
    message('process ', f)
    x = data.table::fread(f, stringsAsFactors = F, data.table = F)
    x$Sample = tools::file_path_sans_ext(basename(f))
    x$Fusion.genes = paste0(x[[1]],"#",x[[2]])
    x <- x %>% dplyr::rename('Hostgene5' = "Gene_1_symbol(5end_fusion_partner)",
                             'Hostgene3' = "Gene_2_symbol(3end_fusion_partner)",
                             'Predicted.effect' = 'Predicted_effect',
                             'Fusion.sequence' = 'Fusion_sequence',
                             'Fusion.description' = 'Fusion_description',
                             'counts' = 'Spanning_unique_reads')
    x <- x[x$counts > min.count,]

    tmp <- data.table::tstrsplit(x$`Fusion_point_for_gene_1(5end_fusion_partner)`, ":", fixed = TRUE)
    x$Chr.5 = tmp[[1]]; x$Position.5 = tmp[[2]]; x$Strand.5 = tmp[[3]]
    tmp <- data.table::tstrsplit(x$`Fusion_point_for_gene_2(3end_fusion_partner)`, ":", fixed = TRUE)
    x$Chr.3 = tmp[[1]]; x$Position.3 = tmp[[2]]; x$Strand.3 = tmp[[3]]

    x = x %>% tidyr::unite(col='chimericRNA', c('Fusion.genes','Chr.5','Position.5','Strand.5','Chr.3','Position.3','Strand.3'),sep="#", remove = F)
    return(x)
  })
  fTable = do.call(rbind, fList)
  fTable$Fusion.seqLen = nchar(fTable$Fusion.sequence)-1
  fTable$Protein.entry.5 = gene_entry$Entry[match(fTable$Hostgene5, gene_entry$Gene.name)]
  fTable$Protein.entry.3 = gene_entry$Entry[match(fTable$Hostgene3, gene_entry$Gene.name)]
  fTable$Fusion.type = ifelse(!is.na(fTable$Protein.entry.5) & fTable$Predicted.effect %in% c('in-frame', 'out-of-frame', 'CDS(truncated)/UTR'), "EE",
                              ifelse(!is.na(fTable$Protein.entry.5), "non-EE", "non-coding"))

  fd = fTable[!duplicated(fTable$chimericRNA),c('chimericRNA', 'Fusion.genes', 'Fusion.description', 'Predicted.effect', 'Fusion.sequence','Fusion.seqLen','Fusion.type',
                                                'Hostgene5', 'Protein.entry.5', 'Chr.5', 'Position.5', 'Strand.5', 'Hostgene3', 'Protein.entry.3', 'Chr.3', 'Position.3', 'Strand.3')]
  fd$Hostgene3 = ifelse(stringr::str_detect(fd$Hostgene3,'^C[0-9]+ORF[0-9]+'),
                        gsub('ORF','orf',fd$Hostgene3), fd$Hostgene3)
  fd$Hostgene5 = ifelse(stringr::str_detect(fd$Hostgene5,'^C[0-9]+ORF[0-9]+'),
                        gsub('ORF','orf',fd$Hostgene5), fd$Hostgene5)

  fd = tibble::remove_rownames(fd)
  fd = tibble::column_to_rownames(fd, var = 'chimericRNA')
  message("import ", nrow(fd), " fusion event.")

  count = reshape2::dcast(fTable, chimericRNA~Sample, value.var ='counts', fill = 0)
  count = tibble::column_to_rownames(tibble::remove_rownames(count), var='chimericRNA')
  count = Matrix::Matrix(Matrix::as.matrix(count), sparse = T)
  count = count[rownames(fd),]

  fd$nSample = Matrix::rowSums(count > 0)
  fd$TotalReads = Matrix::rowSums(count)

  pd = data.frame(nFusion = Matrix::colSums(count > 0), nCount = Matrix::colSums(count))
  obj = new('FusionSet', count = count, sample.data = pd, feature.data = fd, version='0.1.0')

  if(!is.null(info)){
    obj@sample.data[rownames(info), colnames(info)] = info
  }
  return(obj)
}

