#'Rnet_Heatmap - Function to generate a bitmap figure to represent edges in similar networks over time.
#' 
#'Bitmaps can be used to show how edges in a network change over time or other criteria. This function takes an object of the 'rnet.multi.strata' class, specifically the E_matrix, and returns a numerical matrix used to visualize E_matrix. The numerical matrix can be visualized with a call to the 'image' function; The assinged colors are stored in the 'palette' attribute attached to the matrix. This matrix can then be used to plot the bitmap in which each horizontal row in the bitmap represents a unique edge found in the set of networks in rnet.multi.strata object and each vertical column is a network stratum (see the 'Rnet' method for more information). The colors represent the binned values in the E_matrix.
#'
#' @param rnet.list An object of class 'rnet.multi.strata'
#' @param e.cutpoints A vector of numeric values used to cut the edge attribute values in rnet.list@E_matrix. Note binning of negative values is based on absolute value, making the categories symmetric around 0. 
#' @param pos.colors A vector of colors corresponding to the binned positive rnet.list@E_matrix values. Defaults to 4 shades of red.
#' @param neg.colors A vector of colors corresponding to the binned negative rnet.list@E_matrix values. Defaults to 4 shades of green.
#' @param zero.color A single value for coloring edges with value = 0 (typically corresponds to an absent edge).
#' @param NA.color A single value for color invalid edges. Edges will typically be found to be in valid if one or both incident vertices were missing in some strata or insufficient observations were available to estimate an edge in a stratum (see the 'n_threshold' argument in the method 'Rnet').
#' The length of pos.colors and neg.colors must be the same and must be of length 1 greater than e.cutpoints.
#' @export
#' @return The code used to produce the plot.
#' @examples
#' #Example using EC_Rnets_byYear
#' EC_Rnets_byYear <- Rnet(Data = NARMS_EC_DATA, 
#' 						L1 = 0.3, 
#' 						V_set = c('AMP', 'AMC', 'AXO', 'TIO', 'NAL', 
#' 						  'CIP', 'STR', 'GEN', 'COT', 'FIS'), 
#' 						n_threshold = 20,
#'						Stratify = 'Year'
#' 						)
#'
#' EC_Heatmap <- Rnet_Heatmap(EC_Rnets_byYear, e.cutpoints = c(0, 0.05, 0.1, 0.2, 1))
#'
#'par(mar = c(4, 5, 1, 1)+0.1)
#'image(EC_Heatmap, col = attr(EC_Heatmap, 'palette'),
#'	axes = FALSE)

#'axis(1, 
#'	at = seq(0, 1, 1/(dim(EC_Heatmap)[1]-1)),
#'	labels = rownames(EC_Heatmap),
#'	tck = -0.02
#'	)
#'
#'axis(2,
#'	at = seq(0, 1, 1/(dim(EC_Heatmap)[2]-1)),
#'	labels = colnames(EC_Heatmap),
#'	tck = -0.015,
#'	las = 2,
#'	cex.axis = 5/6
#'	)


Rnet_Heatmap <- function(rnet.list,
			e.cutpoints,
			pos.colors = c('#FFBBBB', '#FF8888', '#FF4444', '#FF0000'),
			neg.colors = c('#BBFFBB', '#88FF88', '#44FF44', '#00FF00'),
			zero.color = '#FFFFFF',
			NA.color = '#CCCCCC'
			)
{	
	if(length(pos.colors)!= length(e.cutpoints)-1 | length(neg.colors)!= length(e.cutpoints)-1) stop("vectors pos.colors and neg.colors must be 1 shorter than e.cutpoints")
	names(pos.colors) <- paste('pos', 1:(length(e.cutpoints)-1), sep = '.')
	names(neg.colors) <- paste('neg', 1:(length(e.cutpoints)-1), sep = '.')
	names(zero.color) <- 'zero'
	names(NA.color) <- 'NA'
	palette_set <- c(pos.colors, neg.colors, zero.color, NA.color)

	edge_list <- reshape(as.data.frame(rnet.list@E_matrix),
		direction = 'l',
		idvar = 'Edge',
		ids = rownames(rnet.list@E_matrix),
		varying = list(colnames(rnet.list@E_matrix)),
		v.names = 'Edge_val',
		times = colnames(rnet.list@E_matrix),
		timevar = 'Stratum'
		)

	edge_list$Stratum <- gsub(paste(slot(rnet.list, 'Stratify_by'), '.', sep = ''), '', edge_list$Stratum)
	edge_list$Palette_code <- unlist(sapply(edge_list$Edge_val, function (x) {
		if(is.na(x)) return ('NA')
		switch(
			sign(x) + 2,
			paste('neg', as.integer(cut(abs(x), e.cutpoints)), sep = '.'),	
			'zero',
			paste('pos', as.integer(cut(x, e.cutpoints)), sep = '.'),
			'NA'
			)
		}))

	edge_list$Palette_num <- match(edge_list$Palette_code, names(palette_set))
	color_frame <- reshape(edge_list,
		idvar = 'Edge',
		timevar = 'Stratum',
		drop = c('Edge_val', 'Palette_code'),
		direction =  'w'
		)

	color_mat <- t(as.matrix(color_frame[, -1]))
	dimnames(color_mat) <- list(
		gsub("Palette_num.", '', rownames(color_mat)),
		color_frame$Edge
		)
	attr(color_mat, 'palette') <- palette_set
	
	return(color_mat)
}