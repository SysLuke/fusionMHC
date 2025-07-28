#' calculate differential fusion events
#' @description find differential fusion (DF) events using DESeq2 (based on a model using the negative binomial distribution) and MAST ()
#' @param obj a FusionSet object
#' @param group set the group column
#' @param nameT name of treat group
#' @param nameC name of control group
#' @param at.least.nSample filter the fusion events detected in at least n samples
#' @return results of differential fusion events
#' @export
#' @references
#' Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' McDavid A, Finak G, Yajima M (2023). _MAST: Model-based Analysis of Single Cell Transcriptomics_. doi:10.18129/B9.bioc.
#'
#' @examples
#' \dontrun{
#' df_res = calculate_DF(obj, group='group',nameT='CA',nameC='NAT')
#' head(df_res)
#' }

calculate_DF <- function(obj, group, nameT, nameC, at.least.nSample=1){
  groups = obj@sample.data[[group]]
  stopifnot(c(nameT,nameC) %in% groups)

  counts = Matrix::as.matrix(obj@count)
  samT = colnames(counts)[which(groups == nameT)]
  samC = colnames(counts)[which(groups == nameC)]
  no = setdiff(c(samT, samC), colnames(counts))
  if(length(no) > 0){
    stop("sample(s) not in count matrix: ", paste0(no, collapse = ','), ". Please check them.")
  } else{
    counts = counts[,c(samT,samC)]
  }
  cond = c(rep(nameT,length(samT)), rep(nameC, length(samC)))
  cond = factor(cond, levels = c(nameC,nameT)) # nameC is first/reference
  cond.info = data.frame(row.names = colnames(counts), cond)

  # filter low detection fusion events
  keep = Matrix::rowSums(counts > 0) > at.least.nSample
  counts = counts[keep,]
  message("keep ",nrow(counts), " fusion events.")

  pct.1 = Matrix::rowSums(counts[,samT]>0)
  pct.2 = Matrix::rowSums(counts[,samC]>0)

  # DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(counts+1, colData = cond.info, design=~cond) # pseudo-count of 1 to every entry for calculating geometric mean
  dds2 <- DESeq2::DESeq(dds,quiet = T, fitType='local')
  x <- as.data.frame(DESeq2::results(dds2, contrast=c("cond",nameT,nameC),cooksCutoff=FALSE, independentFiltering=FALSE))
  data.table::setnames(x, old = c('log2FoldChange','P.Value','pvalue','PValue','adj.P.Val','FDR','padj'),
                       new = c('log2FC','Pvalue_deseq2','Pvalue_deseq2','Pvalue_deseq2','Padj_deseq2','Padj_deseq2','Padj_deseq2'), skip_absent=T)
  x = cbind(pct.1=pct.1[rownames(x)],pct.2=pct.2[rownames(x)],x)
  x = tibble::rownames_to_column(x, var='chimericRNA')
  x$baseMean = x$baseMean-1

  # MAST
  latent.vars = cond.info
  latent.vars$wellKey <- rownames(latent.vars)
  fdat <- data.frame(primerid = rownames(counts), row.names = rownames(counts))
  sca <- MAST::FromMatrix(exprsArray = as.matrix(x = counts),
                          check_sanity = FALSE, cData = latent.vars, fData = fdat)
  cdr2 <- Matrix::colSums(counts>0)
  sca@colData$cngeneson <- scale(cdr2)
  zlmfit <- MAST::zlm(formula = ~cond, sca = sca)
  fc = MAST::logFC(zlmfit)
  summaryCond <- MAST::summary(object = zlmfit, doLRT = paste0('cond',nameT))
  summaryDt <- summaryCond$datatable
  y <- summaryDt[which(summaryDt[, "component"] == "H"), c('primerid','Pr(>Chisq)')]
  colnames(y) = c('chimericRNA','Pvalue_mast')
  y$Padj_mast = stats::p.adjust(y$Pvalue,method = 'fdr')
  o = merge(x,y,by='chimericRNA')

  return(o)
}
