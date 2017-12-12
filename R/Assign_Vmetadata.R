#'Assign_Vmetadata - Methods for assigning network vertex metadata
#'
#' This method assigns metadata to vertices in an igraph object based on an existing vertex attibute, typically vertex name, but can be done with other attributes including those assigned with previous calls with this method. 
#'
#' This method also works with all rnet objects (currently class 'rnetBasic', 'rnetStrata', and 'rnet.multi.strata'), and also adds the names of the metadata attributes to the 'V_metadata' slot.
#'
#' @param network The network to which the vertex metadata will be applied.
#' @param V_metadata A dataframe containing the vertex metadata to be assigned. A vertex attribute will be assigned for every column in the frame, except the column used to match V_metadata to existing vertex attribues.
#' @param match.attr The name of the column in V_metadata used to match metadata to vertices. Defaults to the first column of V_metadata.
#' @param vertex.match.attr the name of the vertex attribute used to match metadata. Defaults to 'name' (V(network)$name), which is typically assigned when the network is created with igraph functions.
#' @param reassign A logical argument controling if the function should overwrite the called network argument. Defaults to 'True' for brevity.
#' @import igraph
#' @return An object of the same type as x, with the new vertex attributes assigned by matching 'match.attr' to 'vertex.match.attr'.
#' @rdname Assign_Vmetadata
#' @include Rnet_classes.R
#' @examples 
#' # V_ATTRS' is a data.frame included in the package containing vertex metadata
#' #     regarding antimicrobial class and a color scheme for the vertices. These 
#' #     attributes are useful for plotting and determining modularity.
#' 
#' ABX_LIST <- c('AMP', 'AMC', 'AXO', 'TIO', 'NAL', 'CIP', 'STR', 'GEN', 'COT', 'FIS')
#' 
#' EC08_rnet <- Rnet(NARMS_EC_DATA, 
#'   L1 = 0.25, 
#'   V_set = ABX_LIST, 
#'   Stratify = NARMS_EC_DATA$Year == 2008
#'   )
#' 
#' #Attributes prior to additions
#' vertex_attr_names(EC08_rnet@R)
#' vertex_attr(EC08_rnet@R)
#' 
#' Assign_Vmetadata(EC08_rnet, 
#'                  V_metadata = V_ATTRS,
#'                  match.attr = 'Code',
#'                  vertex.match.attr = 'name'
#'                  )
#' 
#' #NOTE: EC08_rnet does not need to be reassigned for brevity. Returns data.frame of assigned data.
#' #      Reassignment can be performed, if desired. data.frame not returned in such a case.
#' 
#' EC08_withAttrs <- Assign_Vmetadata(EC08_rnet, 
#'                  V_metadata = V_ATTRS,
#'                  match.attr = 'Code',
#'                  vertex.match.attr = 'name'
#'                  )
#' 
#' #Atrributes after edges assigned.
#' vertex_attr(EC08_rnet@R)

#' @export

setGeneric('Assign_Vmetadata',
	function(network, V_metadata, match.attr = NULL, vertex.match.attr = 'name', reassign = T)
	{
		if(is.null(match.attr)) {
			match.attr <- names(V_metadata)[1]
			warning("No column containing vertex names declared, first column of V_metadata is assumed to contain vertex names")
		}								#Assigns first column of metadata table as matching column if none defined and returns warning. 
		if(any(duplicated(V_metadata[match.attr]))) stop("Elements of matching column 'match.attr' in must be unique")
		if(!match.attr%in%names(V_metadata)) stop(paste('Column', match.attr, 'not found in V_metadata'))
										#Returns error if assigned matching column in met

		if(!all(vertex_attr(network, vertex.match.attr)%in%unlist(V_metadata[match.attr]))) stop('Not all network vertices appear in V_metadata$', match.attr, sep = '')
		match.vec <- match(vertex_attr(network, vertex.match.attr), V_metadata[[match.attr]])
		source.env <- parent.frame()

		attr.frame <- data.frame(V = vertex_attr(network, 'name'))

		for(attrib in names(V_metadata)[!names(V_metadata)%in%match.attr]) {
			network <- set_vertex_attr(network, attrib, value = V_metadata[match.vec, attrib])
			attr.frame[[attrib]] <- V_metadata[match.vec, attrib]
		}
#		browser()

		if(reassign) {
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				source.env
				)
			return(attr.frame)
		}
		return(network)
	})

#' @rdname Assign_Vmetadata
#'
setMethod('Assign_Vmetadata',
	signature(network = 'rnetBasic'),
	function(network, V_metadata, match.attr = NULL, vertex.match.attr = 'name', reassign = T)
	{
		network@R <- Assign_Vmetadata(network@R, V_metadata, match.attr, vertex.match.attr, F)
		network@V_metadata <- names(V_metadata)

		if(reassign) {
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(vertex_attr(network@R)))
		}
		return(network)
	})

#' @rdname Assign_Vmetadata
#'
setMethod('Assign_Vmetadata',
	signature(network = 'rnetStrata'),
	function(network, V_metadata, match.attr = NULL, vertex.match.attr = 'name', reassign = T)
	{
		network@R <- Assign_Vmetadata(network@R, V_metadata, match.attr, vertex.match.attr, F)
		network@V_metadata <- names(V_metadata)

		if(reassign) {
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(vertex_attr(network@R)))
		}
		return(network)
	})

#' @rdname Assign_Vmetadata
#'	
setMethod('Assign_Vmetadata',
	signature(network = 'rnetMultiStrata'),
	function(network, V_metadata, match.attr = NULL, vertex.match.attr = 'name', reassign = T)
	{
	
		slot(network, "R_Strata") <- lapply(slot(network, "R_Strata"), Assign_Vmetadata, V_metadata, match.attr, vertex.match.attr, reassign = F)
		
		if(reassign) {
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(vertex_attr(network@R_Strata[[1]]@R)))
		}
		return(network)
	})