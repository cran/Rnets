 #'.Gen_R - internal methods for generating Rnets.
#'
#'Internal for the "Rnet" methods. .Gen_R should not be called directly.
#' @param rnet.obj an Rnet S4 object created by Rnet()
#' @importFrom igraph graph_from_adjacency_matrix
#' @import glasso
#' @importFrom stats reshape
#' @rdname dot-Gen_R
#' 
setGeneric('.Gen_R', function(rnet.obj){

	rnet.obj@Data <- rnet.obj@RawData[rnet.obj@V_set_orig]
	rnet.obj@L1 <- rnet.obj@L1_orig

	if(!all(rnet.obj@V_set_orig%in%names(rnet.obj@Data))) stop("One or more vertex names declared in 'V_set' do not match data column names.")

	rnet.obj@V_set <- rnet.obj@V_set_orig[lapply(rnet.obj@Data, function(x) sum(!is.na(x))) > rnet.obj@n_threshold]
	V_sorted <- sort(rnet.obj@V_set)
	rnet.obj@Data <- rnet.obj @Data[V_sorted]

	if(!all(apply(rnet.obj@Data, c(1, 2), is.numeric))) stop("Non-numeric results provided in data for correlation matrix")

	rnet.obj@Sigma <- Clean_Sigma(rnet.obj@Data, rnet.obj@cor_method, rnet.obj@cor_pairing)

	rnet.obj@n <- t(!is.na(rnet.obj@Data))%*%!is.na(rnet.obj@Data)
	n_list <- Sq2L(rnet.obj@n, c("V1", "V2", "n"), drop.values = NULL)
	rnet.obj@Zeros <- list(
		Forced = if(dim(rnet.obj@Forced_zeros)[1]!=0) rnet.obj@Forced_zeros else matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c('V1', 'V2'))),
		Invalid = data.matrix(subset(n_list, select = c(V1, V2), n < rnet.obj@n_threshold))
		)
	zero_pairs <- unique(rbind(rnet.obj@Zeros$Forced, rnet.obj@Zeros$Invalid))
	zero_indices <- data.matrix(data.frame(V1 = match(zero_pairs[,1], V_sorted), V2 = match(zero_pairs[,2], V_sorted)))

	if(dim(zero_pairs)[1]==0) 	glasso_result <- glasso(rnet.obj@Sigma, rnet.obj@L1) else glasso_result <- glasso(rnet.obj@Sigma, rnet.obj@L1, zero = zero_indices)

	rnet.obj@Theta <- glasso_result$wi
	dimnames(rnet.obj@Theta) <- list(V_sorted, V_sorted)

	rnet.obj@Omega <- Estimate_pCorrs(rnet.obj@Theta)

	rnet.obj@A <- rnet.obj@Omega!=0

	rnet.obj@R <-  graph_from_adjacency_matrix(rnet.obj@A, mode = 'undirected')
	
	E(rnet.obj@R)$omega <- merge(x = as.data.frame(
	                               as_edgelist(rnet.obj@R), 
	                               stringsAsFactors = F
	                               ), 
	                             y = Sq2L(
	                               rnet.obj@Omega, 
	                               c("V1", "V2", "omega"), 
	                               c(T, F, F),
	                               0
	                               ),
	                             by = c('V1', 'V2')
	                             )$omega
	
	rnet.obj@E_metadata <- edge_attr_names(rnet.obj@R)
	
  if(is.null(rnet.obj@Layout_master)) rnet.obj@Layout <- layout_with_fr(rnet.obj@R)
	
	rnet.obj@V_omitted <- rnet.obj@V_set_orig[!rnet.obj@V_set_orig%in%rnet.obj@V_set]

	return(rnet.obj)
})

#'Internal for the "Rnet" methods. .Gen_R should not be called directly.

#' @rdname dot-Gen_R

#' 
setMethod('.Gen_R',
	'rnetStrata',
	function (rnet.obj) {
		rnet.obj@Data <- rnet.obj@RawData[eval(rnet.obj@Strata_def, rnet.obj@RawData),rnet.obj@V_set_orig]
		rnet.obj@L1 <- rnet.obj@L1_orig

		if(!all(rnet.obj@V_set_orig%in%names(rnet.obj@Data))) stop("One or more vertex names declared in 'V_set' do not match data column names.")

		rnet.obj@V_set <- rnet.obj@V_set_orig[lapply(rnet.obj@Data, function(x) sum(!is.na(x))) > rnet.obj@n_threshold]
		V_sorted <- sort(rnet.obj@V_set)
		rnet.obj@Data <- rnet.obj@Data[V_sorted]

		if(!all(apply(rnet.obj@Data, c(1, 2), is.numeric))) stop("Non-numeric results provided in data for correlation matrix")

		rnet.obj@Sigma <- Clean_Sigma(rnet.obj@Data, rnet.obj@cor_method, rnet.obj@cor_pairing)

		rnet.obj@n <- t(!is.na(rnet.obj@Data))%*%!is.na(rnet.obj@Data)
		n_list <- Sq2L(rnet.obj@n, c("V1", "V2", "n"), drop.values = NULL)

		rnet.obj@Zeros <- list(
			Forced = if(dim(rnet.obj@Forced_zeros)[1]!=0) rnet.obj@Forced_zeros else matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c('V1', 'V2'))),
			Invalid = data.matrix(subset(n_list, select = c('V1', 'V2'), n_list$n < rnet.obj@n_threshold))
			)
		zero_pairs <- unique(rbind(rnet.obj@Zeros$Forced, rnet.obj@Zeros$Invalid))
		zero_indices <- data.matrix(data.frame(V1 = match(zero_pairs[,1], V_sorted), V2 = match(zero_pairs[,2], V_sorted)))


		if(dim(zero_pairs)[1]==0) 	glasso_result <- glasso(rnet.obj@Sigma, rnet.obj@L1) else glasso_result <- glasso(rnet.obj@Sigma, rnet.obj@L1, zero = zero_indices)

		rnet.obj@Theta <- glasso_result$wi
		dimnames(rnet.obj@Theta) <- list(V_sorted, V_sorted)
	
		rnet.obj@Omega <- Estimate_pCorrs(rnet.obj@Theta)
		
		rnet.obj@A <- rnet.obj@Omega!=0
		
		rnet.obj@R <- graph.adjacency(rnet.obj@A, mode = 'undirected')

		E(rnet.obj@R)$omega <- merge(
		  x = as.data.frame(
		    as_edgelist(rnet.obj@R), 
		    stringsAsFactors = F
		  ), 
		  y = Sq2L(
		    rnet.obj@Omega, 
		    c("V1", "V2", "omega"), 
		    c(T, F, F),
		    0
		  ),
		  by = c('V1', 'V2')
		)$omega
		rnet.obj@A <- rnet.obj@Omega!=0		
		if(is.null(rnet.obj@Layout_master)) rnet.obj@Layout <- layout_with_fr(rnet.obj@R)
		
		rnet.obj@V_omitted <- rnet.obj@V_set_orig[!rnet.obj@V_set_orig%in%rnet.obj@V_set]

		return(rnet.obj)
	})