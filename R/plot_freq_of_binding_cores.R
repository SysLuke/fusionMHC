#' plot top frequent MHC binding cores
#' @description barplot shows the most frequent MHC binding cores
#' @param obj a FusionSet object
#' @param topN show top N binding cores
#' @param group set the group name
#' @param hla.type set the HLA type
#' @param point.size set point.size
#' @param bar.color set bar.color
#' @param bar.width set bar.width
#' @return a ggplot2 object
#'
#' @import ggplot2 dplyr
#' @importFrom utils tail
#' @export
#' @examples
#' \dontrun{
#' top_cores = plot_most_freq_binding_cores(obj, topN=10, group="group")
#' }

plot_most_freq_binding_cores <- function(obj, topN = 10, group = NULL, hla.type = c('HLA-A0201','HLA-A1101'),
                                        point.size = 5, bar.color = 'gray', bar.width = 0.1){
  pt = obj@peptide
  if(nrow(pt) == 0) stop('no records. please check it')
  if(!all(c('HLA-A0201:aaseq',"HLA-A1101:aaseq") %in% names(obj@peptide))){
    stop("no prediction of binding. please run predict_MHC_affinity first.")
  }

  pl = list()
  for(ht in hla.type){
    cn = paste0(ht,':aaseq')
    tmp = pt[which(!is.na(pt[cn])),]

    if(is.null(group)){
      cnt = sort(table(tmp[cn]))
      cnt = tail(cnt, min(topN,length(cnt)))
      cntdf = as.data.frame(cnt)
      colnames(cntdf)[1] = 'binding cores'
      out = tmp[tmp[[cn]] %in% names(cnt),]
      g = ggplot(cntdf, mapping = aes(x=.data[['Freq']], y=.data[['binding cores']])) +
        geom_col(width = bar.width, fill = bar.color) +
        geom_point(aes(color = .data[['Freq']]),size=point.size) +
        ggpubr::theme_classic2(base_size = 14) +
        guides(color = 'none') + labs(title = ht) + ylab('binding cores') +
        scale_color_distiller(palette = "Spectral", direction = -1) +
        scale_x_continuous(limits = c(NA,max(cnt)+0.1*max(cnt)))
    } else{
      stopifnot(group %in% names(obj@sample.data))

      pd = obj@sample.data
      counts = obj@count

      sg = split(pd, f=pd[[group]])

      pl2 = list()
      out = data.frame()
      for(i in names(sg)){
        x = sg[[i]]
        counts_i = counts[,rownames(x)]
        cnt_i = Matrix::rowSums(counts_i)
        fts = names(cnt_i)[cnt_i>0]
        tmp = pt[match(fts,pt$chimericRNA),]

        cnt = sort(table(tmp[cn]))
        cnt = tail(cnt, min(topN,length(cnt)))
        cntdf = as.data.frame(cnt)
        colnames(cntdf)[1] = 'binding cores'
        cntdf$group = i

        tmp = tmp[tmp[[cn]] %in% names(cnt),]
        tmp$group = i
        out = rbind(out,tmp)

        g = ggplot(cntdf, mapping = aes(x=.data[['Freq']], y=.data[['binding cores']])) +
          geom_col(width = bar.width, fill = bar.color) +
          geom_point(aes(color = .data[['Freq']]),size=point.size) +
          ggpubr::theme_classic2(base_size = 14) +
          guides(color = 'none') + labs(title = ht, subtitle = paste0('group:',i)) + ylab('binding cores') +
          scale_color_distiller(palette = "Spectral", direction = -1) +
          scale_x_continuous(limits = c(NA,max(cnt)+0.1*max(cnt)))
        pl2[[i]] = g
      }
      g = cowplot::plot_grid(plotlist = pl2)
    }
    pl[[ht]] = g
  }
  g = cowplot::plot_grid(plotlist = pl)
  print(g)
  return(out)
}
