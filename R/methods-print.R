#'Print methods of rnet objects.
#'
#'print methods for classes 'rnetBasic', 'rnetStrata', 'rnet.multi.strata'.
#' @param x An rnet object
#' @rdname print
#' @export
#' 
setMethod(f = 'print',
          signature(x = 'rnetBasic'),
          function(x) {
            cat(  '\nVertex set:',	paste(x@V_set, collapse = ' '),
                  '\n',
                  '\n  Edge Set:',
                  '\n')
            edge.frame <- data.frame(V1 = as_edgelist(x@R)[,1],
                                     V2 = as_edgelist(x@R)[,2],
                                     Omega = round(E(x@R)$omega, 3)	
            )
            print(edge.frame)
          })

#' @rdname print
#' 
setMethod(f = 'print',
          signature(x = 'rnetStrata'),
          function(x) {
            cat(	'\nRnet conditioned on', as.character(x@Strata_def), '\n')
            print(as(x, 'rnetBasic'))
          })

#' @rdname print
#' 
setMethod(f = 'print',
          signature(x = 'rnetMultiStrata'),
          function(x) {
            cat(	"\nRnet stratified by ", x@Stratify_by, " (", length(x@R_Strata), " Strata)\n", sep = '')
            cat(  '\nVertex set:',	paste(unique(unlist(lapply(x@R_Strata, function(x) x@V_set ))), collapse = ' '), "\n")
            cat(  '\nEdge set:\n')
            print(.Assemble_Edge_Matrix(x, 'omega'))
            
          })
