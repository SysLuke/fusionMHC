# class definition for class FusionSet
#' @slot count reads count
#' @slot sample.data data frame storing the sample information
#' @slot feature.data data frame storing the meta information of fusion events
#' @slot peptide prediction of binding to MHC class I molecules
#' @slot PTC prediction of PTC caused by NMD
#' @slot version package version
#' @importFrom methods setClass setMethod
FusionSet <- setClass(Class = 'FusionSet',
                      slots = c(count = 'dgCMatrix',
                               occurrence = 'list',
                               sample.data = 'data.frame',
                               feature.data = 'data.frame',
                               peptide = 'data.frame',
                               PTC = 'data.frame',
                               version = 'character'))

setMethod(
  f = 'show', signature = signature(object='FusionSet'),
  definition = function(object) {
    cat('An object of class ',
        class(object), ': ',
        nrow(object@count), ' fusion event across ',
        ncol(object@count), ' samples.\n',
        sep = '')
    invisible(F)
  }
)
