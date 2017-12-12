#'Assign_Emetadata - Methods for assigning network edge metadata
#'
#' This method assigns metadata to edges in an igraph object based on a continuous edge attribute already present in object, typically called with E(graph)$attr. This method is designed to be a convenient alternative to repeated calls to 'E(graph)$attr <-  value' or 'set_edge_attr' to apply multiple attributes that vary on a single edge attribute. The attribute used to assign the other metadata attributes is typically assumed to be the value used to define the edge to begin with, so when the attribute is 0, it is typically assumed the edge doesn't exist and is omitted from the edge set.
#'
#' This method also works with all rnet objects (currently class 'rnet.basic', 'rnet.strata', and 'rnet.multi.strata'), and also adds the names of the metadata attributes to the 'E_metadata' slot.

#' @param network The network to which the edge metadata will be applied.
#' @param E_metadata A dataframe containing one column for each attribute to be assigned to the network edges.
#' @param match.attr The continuous attribute assigned to the network's edges (typically by E(network)$attr or set_edge_attr(network, attr)) that will be used to categorize edges. Values in match.attr will be binned using cut(value, E_ctpts), and matched to the rows in E_metadata, i.e. edges binned into the first category will be assigned attribute values from the first row in the data set, those binned into the second category will be assigned the second row values, etc.
#' @param e.cutpoints A vector of values used to define how the values of match.attr will be binned. See 'cut' for more information.
#' @param sign.color By default, a two-element vector containing colors for edges representing positive and negative attribute values(black and red, respectively). This behavior will be overriden if this argument is set to FALSE, NA, or NULL. Note, edge_attr('color') will ALWAYS be assinged to network edges when this method is called.
#' @param attr_abs_val A logical argument determining if the absolute value of match.attr is used when binning values.
#' @param reassign A logical argument controling if the function should overwrite the called network argument. Defaults to 'True' for brevity.
#' @import igraph
#' @include Rnet_classes.R
#' @return An object of the same type as x, with the new edge attributes assigned by binning the 'match.attr' in the igraph and assigning the matching rows in 'e.metadata'.

#' @export
#' @examples
#'
#' #'E_ATTRS' is a data.frame included in the package containing line type and weight  
#' #     for plotting edges ('lty' & 'width', respectively) that represent partial
#' #     correlation magnitude and sign held in attribute 'omega'. 
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
#' edge_attr_names(EC08_rnet@R)
#' edge_attr(EC08_rnet@R)
#' 
#' OMEGA_CUTS <- c(0, 0.05, 0.10, 0.20, 1) #Cutpoints to sort abs(omega) into 4 bins
#' 
#' Assign_Emetadata(EC08_rnet, 
#'                  E_metadata = E_ATTRS,
#'                  match.attr = 'omega',
#'                  e.cutpoints = OMEGA_CUTS
#'                  )
#' 
#' #NOTE: EC08_rnet does not need to be reassigned for brevity. Returns data.frame of assigned data.
#' #      Reassignment can be performed, if desired. data.frame not returned in such a case.
#' 
#' EC08_withAttrs <- Assign_Emetadata(EC08_rnet, 
#'                  E_metadata = E_ATTRS,
#'                  match.attr = 'omega',
#'                  e.cutpoints = OMEGA_CUTS
#'                  )
#' 
#' #Atrributes after edges assigned.
#' edge_attr(EC08_rnet@R)
#' 
#' #NOTE: Edge color assigned as per default behavior.
#' @rdname Assign_Emetadata
setGeneric('Assign_Emetadata',
	function(network, E_metadata, match.attr, e.cutpoints = NULL, sign.color = c('black', 'red'), attr_abs_val = T, reassign = T)
	{
		if(!match.attr%in%edge_attr_names(network)) stop(paste("Edge attribute '", match.attr, "' not found in network", sep = ''))
		if(!is.numeric(edge_attr(network, match.attr))) stop("match.attr is not numeric")
		if(attr_abs_val) match.attr_vals <- abs(edge_attr(network, match.attr)) else match.attr_vals <- edge_attr(network, match.attr)
		if(is.null(e.cutpoints)) e.cutpoints <- seq(min(0, min(match.attr_vals)), max(match.attr_vals), max(match.attr_vals)/dim(E_metadata)[1])
		E_cat <- cut(match.attr_vals, e.cutpoints)

		for(attr in names(E_metadata)) network <- set_edge_attr(network, attr, value = E_metadata[E_cat,attr])
		if(!'color'%in%names(E_metadata)) if(sign.color[1] == FALSE|is.null(sign.color)[1]|is.na(sign.color)[1]) E(network)$color <- 'black' else {
			E(network)$color <- sign.color[1]; E(network)$color[sign(edge_attr(network, match.attr))==-1] <- sign.color[2]
			}

		if(reassign){
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(edge_attr(network@R)))
		}

		
		return(network)
	})

#' @rdname Assign_Emetadata
#' 
setMethod('Assign_Emetadata',
	signature(network = 'rnetBasic'),
	function(network, E_metadata, match.attr, e.cutpoints = NULL, sign.color = c('black', 'red'), attr_abs_val = T, reassign = T) 
	{
		network@R <- Assign_Emetadata(network@R, E_metadata, match.attr, e.cutpoints, sign.color, attr_abs_val, F)
		network@E_metadata <- names(E_metadata)

		if(reassign){
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(edge_attr(network@R)))
		}
		return(network) 

	})

#' @rdname Assign_Emetadata
#' 
setMethod('Assign_Emetadata',
	signature(network = 'rnetMultiStrata'),
	function(network, E_metadata, match.attr, e.cutpoints = NULL, sign.color = c('black', 'red'), attr_abs_val = T, reassign = TRUE) 
	{	
		slot(network, "R_Strata") <- lapply(slot(network, "R_Strata"), Assign_Emetadata, E_metadata, match.attr, e.cutpoints, sign.color, attr_abs_val, FALSE)

		if(reassign){
			assign(as.character(as.list(sys.call())[[2]]),
				network,
				parent.frame()
				)
			return(as.data.frame(edge_attr(network@R_Strata[[1]]@R)))
		}
		return(network) 
	})