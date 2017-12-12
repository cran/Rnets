#' Plot methods for R-nets
#'
#' A plot method for R-nets, and incorporates vertex and edge metadata and layout, if assigns. Only vertex and edge metadata with names that match igraph decoration options (without 'vertex.' or 'edge.' appended to the attribute name; see plot.igraph). Layout is pulled from 'Layout_master' in the rnet object, if it exists. The layout frame can contain 3 columns, with the first column used to match the coordinates in the next two columns to graph vertices OR can contain 2 columns with the same number of vertices in the graph.
#' @param x an rnet object of class 'rnet.basic'
#' @rdname plot-Rnet.basic
#' @aliases plot
#' @export

setMethod('plot',
	signature(x = 'rnetBasic'),
	function (x) {
		VERT.PARAMS <- c('size','size2','color','frame.color','shape','label','label.family','label.font','label.cex','label.dist','label.degree','label.color')
		EDGE.PARAMS <- c('color','width','lty','label','label.family','label.font','label.cex','label.color','label.x','label.y','curved')

		vert.attrib.lines <- character(0)
		vert.attribs <- intersect(x@V_metadata, VERT.PARAMS)
		if(length(vert.attribs) > 0) for(attrib.name in vert.attribs) vert.attrib.lines <- c(vert.attrib.lines, paste('vertex.', attrib.name, '= V(x@R)$', attrib.name, sep = ''))

		edge.attrib.lines <- character(0)
		edge.attribs <- intersect(x@E_metadata, EDGE.PARAMS)
		if(length(edge.attribs) > 0) for(attrib.name in edge.attribs) edge.attrib.lines <- c(edge.attrib.lines, paste('edge.', attrib.name, '= E(x@R)$', attrib.name, sep = ''))

		x@Layout <- Rnets::.Assign_Layout_Matrix(x)
		attrib.lines <- c(vert.attrib.lines, edge.attrib.lines)

		if(length(attrib.lines)==0) plot.call <- 'plot.igraph(x@R, layout = x@Layout)' else plot.call <- paste('plot.igraph(x@R', paste(vert.attrib.lines, collapse = ','), 'layout = x@Layout)', sep = ',')

		#rnetBasic()

		eval(parse(text = plot.call))
		#plot.call <- gsub('x', as.list(sys.call(1))[[2]], plot.call)
		#return(plot.call)
	})

#' Plot method for rnetStrata class
#' @param x an rnet object of class 'rnetStrata'
#' @rdname plot-rnetStrata
#' 
setMethod('plot',
	signature(x = 'rnetStrata'),
	function(x) {
		VERT.PARAMS <- c('size','size2','color','frame.color','shape','label','label.family','label.font','label.cex','label.dist','label.degree','label.color')
		EDGE.PARAMS <- c('color','width','lty','label','label.family','label.font','label.cex','label.color','label.x','label.y','curved')

		vert.attrib.lines <- character(0)
		vert.attribs <- intersect(x@V_metadata, VERT.PARAMS)
		if(length(vert.attribs) > 0) for(attrib.name in vert.attribs) vert.attrib.lines <- c(vert.attrib.lines, paste('vertex.', attrib.name, '= V(x@R)$', attrib.name, sep = ''))

		edge.attrib.lines <- character(0)
		edge.attribs <- intersect(x@E_metadata, EDGE.PARAMS)
		if(length(edge.attribs) > 0) for(attrib.name in edge.attribs) edge.attrib.lines <- c(edge.attrib.lines, paste('edge.', attrib.name, '= E(x@R)$', attrib.name, sep = ''))

		x@Layout <- Rnets::.Assign_Layout_Matrix(x)
		attrib.lines <- c(vert.attrib.lines, edge.attrib.lines)
		
		if(length(attrib.lines)==0) plot.call <- 'plot.igraph(x@R, layout = x@Layout)' else plot.call <- paste('plot.igraph(x@R', paste(vert.attrib.lines, collapse = ','), 'layout = x@Layout)', sep = ',')
		
		#rnetBasic()
		
		eval(parse(text = plot.call))
		#plot.call <- gsub('x', as.list(sys.call(1))[[2]], plot.call)
		#return(plot.call)
	})

#' Hidden function for assigning layout matrix
#' @param x an Rnet object
#' @rdname dot-Assign_Layout_Matrix
.Assign_Layout_Matrix <- function(x){ 
  if(is.null(x@Layout_master)) return(layout_with_fr(x@R))
	
  if(dim(x@Layout_master)[2] == 3) {
    coord_ref.vec <- x@Layout_master[,1]
    coord_x.vec <- x@Layout_master[,2]
    coord_y.vec <- x@Layout_master[,3]
  } else if (dim(x@Layout_master)[2] == 2){
    coord_ref.vec <- NULL
    coord_x.vec <- x@Layout_master[,1]
    coord_y.vec <- x@Layout_master[,2]
  } else stop('Layout_master invalid must have 2 or 3 column')

  if(!is.numeric(coord_x.vec)|!is.numeric(coord_y.vec)) stop('Layout_master not valid: Coordinate columns must be numeric')
  
  if(is.factor(coord_ref.vec)) coord_ref.vec <- as.character(coord_ref.vec)
  if(is.null(coord_ref.vec)) {
    if(length(x@V_set)!= length(coord_x.vec)) stop('Layout_master not valid: Number of rows in layout frame must match number of vertices in graph OR a column for vertex matching must be provided.')
    layout.mat <- cbind(coord_x.vec, coord_y.vec)

    dimnames(layout.mat) <- list(V(x@R)$name,c('x', 'y'))
    } else {
      if(!all(V(x@R)$name%in%coord_ref.vec)) stop('Layout_master not valid: Vertex names are missing from the first column.')
      seq.vec <- match(V(x@R)$name,coord_ref.vec)
	    layout.mat <- cbind(coord_x.vec[seq.vec], coord_y.vec[seq.vec])
	    dimnames(layout.mat) <- list(V(x@R)$name[seq.vec],c('x', 'y'))
    }   
	return(layout.mat)
	}