#' plot basic summary
#' @description plot basic summary
#' @param obj a FusionSet object
#' @param features colnames in obj@sample.data
#' @param group specify group name
#' @param group.colors set group colors
#' @param plot.type boxplot or barplot
#' @param topN show top N samples
#' @param bar.width set bar.width
#' @param ncol number of columns to combine multiple feature plots to
#' @return a ggplot2 object
#' @import ggplot2 ggpubr
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' obj = plot_sample_summary(obj,features = c('nFusion','nCount'), group='group')
#' }

plot_sample_summary <- function(obj,
                                features = c('nFusion','nCount'),
                                group = NULL,
                                group.colors = NULL,
                                plot.type = "barplot",
                                topN = 50,
                                bar.width = 0.5,
                                ncol = 2){
  pd = tibble::rownames_to_column(obj@sample.data, var='sample')
  stopifnot(features %in% names(pd))

  if(plot.type == 'barplot'){
    pl = list()
    for(f in features){
      tmp = pd[order(pd[[f]],decreasing = T),]
      tmp = tmp[1:min(nrow(tmp), topN),]
      tmp$sample = factor(tmp$sample, levels = tmp$sample)

       g = ggplot(data = tmp, mapping = aes(x=.data[['sample']], y=.data[[f]], fill = .data[[f]])) +
         ggpubr::theme_classic2(base_size = 14) +
        geom_bar(width = bar.width, stat = 'identity') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = unit(c(5,5,5,10),'pt'),
              legend.key.width = unit(10, 'pt')) + labs(x=NULL, y=NULL, title = f) +
        guides(fill = guide_colorbar(title = NULL)) +
        scale_fill_distiller(palette = "Spectral", direction = -1) +
         guides(fill = 'none')

       if(!is.null(group)){
         tmp$.tmp = tmp[[group]]
         g = ggplot(data = tmp, mapping = aes(x=.data[['sample']], y=.data[[f]], fill = .data[[f]])) +
           ggpubr::theme_classic2(base_size = 14) +
           geom_bar(width = bar.width, stat = 'identity') +
           theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = unit(c(5,5,5,10),'pt'),
                 legend.key.width = unit(10, 'pt')) + labs(x=NULL, y=NULL, title = f) +
           guides(fill = guide_colorbar(title = NULL)) +
           scale_fill_distiller(palette = "Spectral", direction = -1) +
           facet_wrap(~.tmp, scales = 'free') +
           theme(strip.background = element_rect(colour = 'transparent')) +
           guides(fill = 'none')
         tmp$.tmp = NULL
       }
       pl[[f]] = g
    }
    g = cowplot::plot_grid(plotlist = pl, ncol=ncol)
  } else if(plot.type == 'boxplot'){
    if(is.null(group)) stop('it requires group for boxplot')

    pl = list()
    for(f in features){
      if(is.null(group.colors)){
        ngroup = length(unique(pd[[group]]))
        group.colors = setNames(grDevices::rainbow(ngroup), unique(tmp[[group]]))
      }
      g = ggplot(data = pd, mapping = aes(x = .data[[group]], y = .data[[f]], fill = .data[[group]])) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        labs(x=NULL, y=NULL, title = f) +
        scale_fill_manual(values = group.colors) +
        theme_classic2(base_size = 14) + guides(fill = 'none')
      pl[[f]] = g
    }
    g = cowplot::plot_grid(plotlist = pl, ncol=ncol)
  }
  return(g)
}
