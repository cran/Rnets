#'Summary - rnetBasic
#'
#'Gives more information than 'print'.
#' @param object an rnet object of class 'rnetBasic'
#' @param ... Additional arguments passed to 'summary' method
#' @rdname summary-rnetBasic
#' @aliases summary
#' @export

setMethod(f = 'summary',
          signature(object = 'rnetBasic'),
          function(object) {
            cat(	'\nBasic R-net',
                 '\n',
                 '\n Sample:', dim(object@Data)[1], 'isolates,', length(object@V_set), 'vertices ')
            if(length(object@V_omitted)==1) cat('(', length(object@V_omitted), ' vertex omitted)', sep = '') 		
            if(length(object@V_omitted) > 1) cat('(', length(object@V_omitted), ' vertices omitted)', sep = '') 		
            cat(	'\n',
                 '\n     L1:', object@L1,
                 '\n  Edges:', ecount(object@R),
                 '\nDensity:', round(100*edge_density(object@R), 1),'%',
                 '\n')
            print(object)
            if(length(object@V_omitted)>0) cat('\nOmitted Vertices:', paste(object@V_omitted, collapse = ', '), '\n')
            cat('\n')
          })	
#'Summary - rnetStrata
#'
#'Gives more information than 'print'.
#' @param object an rnet object of class 'rnetStrata'

#' @rdname summary-rnetStrata

setMethod(f = 'summary',
          signature(object = 'rnetStrata'),
          function(object) {
            cat(  '\nStratfied R-net (single level)',
                  '\n',
                  '\n Sample:', dim(object@Data)[1], 'isolates,', length(object@V_set), 'vertices ')
            if(length(object@V_omitted)==1) cat('(', length(object@V_omitted), ' vertex omitted)', sep = '') 		
            if(length(object@V_omitted) > 1) cat('(', length(object@V_omitted), ' vertices omitted)', sep = '') 		
            cat(	'\n',
                 '\n     L1:', object@L1,
                 '\n  Edges:', ecount(object@R),
                 '\nDensity:', round(100*edge_density(object@R), 1),'%',
                 '\n')
            print(object)
            if(length(object@V_omitted)>0) cat('\nOmitted Vertices:', paste(object@V_omitted, collapse = ', '), '\n')
            cat('\n')
          })