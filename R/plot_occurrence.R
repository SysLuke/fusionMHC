#' heatmap shows occurrence for each HLA type
#' @description heatmap shows occurrence for each HLA type
#' @param obj a FusionSet object
#' @param hla.type set the HLA type
#' @param anno.sample data frame to annotate the samples
#' @param anno.colors a list of colors
#' @param col set 3 colors if hla.type is "both". otherwise 2 colors.
#' @param fusion.to.show set fusion names to show (chimericRNA)
#' @param top.to.show show n fusion events
#' @param show_column_names TRUE or FALSE
#' @param show_row_names TRUE or FALSE
#' @param ... arguments passed to ComplexHeatmap::Heatmap
#' @return a FusionSet object
#'
#' @export
#' @examples
#' \dontrun{
#' obj = plot_occurrence(obj, hla.type = 'both', col = c("blue",'green','yellow'),
#'                       anno.sample = obj@sample.data,
#'                       anno.colors = list(group = c("CA"="red","NAT"='green')),
#'                       top.to.show = 10)
#' }

plot_occurrence <- function(obj, hla.type = c('HLA-A0201','HLA-A1101','both'), anno.sample = NULL,
                            anno.colors = NULL, col = c("blue",'green','yellow'), fusion.to.show = NULL, top.to.show = 10,
                            show_column_names = T, show_row_names = T, ...){
  ha <- NULL
  if(!is.null(anno.sample)){
    ha <- ComplexHeatmap::HeatmapAnnotation(df = anno.sample, col = anno.colors, show_annotation_name = TRUE,
                                            annotation_name_gp = grid::gpar(fontsize = 7))
  }
  for(hla in hla.type){
    if(hla %in% c('HLA-A0201','HLA-A1101')){
      mat = obj@occurrence[[hla]]
      mat = Matrix::as.matrix(mat)
      cnt = Matrix::rowSums(mat>0)
      if(!is.null(fusion.to.show)){
        names.to.show <- fusion.to.show
      } else{
        names.to.show = names(sort(cnt, decreasing = T)[1:top.to.show])
      }
      cnt = cnt[names.to.show]
      row_ha = ComplexHeatmap::rowAnnotation(Freq = ComplexHeatmap::anno_barplot(cnt, border = F, add_numbers = T))

      mat = mat[names.to.show,]
      dt = ComplexHeatmap::Heatmap(mat, top_annotation = ha,
                              show_heatmap_legend = F,
                              column_title = paste0("fusion events binding to ", hla),
                              column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                              col = col,
                              row_names_side = "left",
                              right_annotation = row_ha,
                              show_column_names = show_column_names,
                              show_row_names = show_row_names, ...)
      ComplexHeatmap::draw(dt)
    } else if(hla == 'both'){
    matA = as.data.frame(obj@occurrence$`HLA-A0201`)
    matB = as.data.frame(obj@occurrence$`HLA-A1101`)

    ofs = unique(rownames(matA),rownames(matB))
    oA = matA[ofs,]
    oB = matB[ofs,colnames(matA)]

    rownames(oA) = ofs
    rownames(oB) = ofs

    oA[is.na(oA)] = 0
    oB[is.na(oB)] = 0

    o = oA+oB
    o[oA == 1 & o ==1] = "HLA-A0201"
    o[oB == 1 & o ==1] = "HLA-A1101"
    o[o ==2] = "both"
    o[o ==0] = NA
    cnt = Matrix::rowSums(!is.na(o))

    if(!is.null(fusion.to.show)){
      names.to.show <- fusion.to.show
    } else{
      names.to.show = names(sort(cnt, decreasing = T)[1:top.to.show])
    }
    o = o[names.to.show,]
    cnt = cnt[names.to.show]
    row_ha = ComplexHeatmap::rowAnnotation(Freq = ComplexHeatmap::anno_barplot(cnt, border = F, add_numbers = T))
    dt = ComplexHeatmap::Heatmap(o, col=setNames(col, c('HLA-A0201','HLA-A1101','both')), row_names_side = "left",
                            top_annotation = ha, name = " ", right_annotation = row_ha,
                            show_column_names = show_column_names, column_title = paste0('occurrence events'),
                            show_row_names = show_row_names,...)

    ComplexHeatmap::draw(dt)
    }
  }
}
