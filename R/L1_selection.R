#' L1 Selection for Rnets
#'
#' An implementation of the Stability Approach to Regularization Selection (StARS) method for L1 penalty selection for use with Rnets method.
#' @param Data The dataset containing the MICs
#' @param L1_set The set of candidate L1 penalties to be evaluated for creating a sparse precision matrix. Must be non-negative.
#' @param B The number of subsamples to evaluate network stability. Defaults to 100 subsamples.
#' @param n_b The size of the subsample to be drawn from the data. If 0 < n_b < 1, this is interpreted as a proportion of the data set; if n_b > 1, it is interpreted as a set sample size. Defaults to 50\% of sample size.
#' @param V_set A character vector corresponding to the names of the antibiotics to include in the Rnet. Defaults to an empty list, in which case all columns in 'MIC_data' will be included in the Rnet
#' @param n_threshold The minimum number of observations required for an an estimated correlation to be valid. Defaults to 0, in which case any number of observations will be sufficient to estimate a valid correlation
#' @param cor_method The method used to estimate the correlation matrix. Must be 'pearson', 'spearman', or 'kendall'. Partial matches allowed. Defaults to 'spearman'.
#' @param cor_pairing The method used to determine how NAs are handled when determining which pairs are used to estimate correlations. See 'cor' function documentation for additional information.
#' @param Stratify The rule for stratifying the data, if desired. 
#' @param Forced_zeros Edges to be omitted from the Rnet.
#' @return A vector of D statistics, corresponding the tested L1 values.
#' @import igraph 
#' @import data.table
#' @export
#' @examples 
#' \donttest{
#'  EC_all_L1Selection <- L1Selection(
#'                              Data = NARMS_EC_DATA, 
#'                              L1_set = seq(0.05, 0.50, 0.05),
#'                              n_b = 1500,
#'                              V_set = ABX_LIST
#'                              )

#' round(EC_all_L1Selection@StARS_D, 4)
#' }

setGeneric('L1Selection',
	function(
		Data,
		L1_set,
		B = 100,
		n_b = 0.5,
		V_set = NULL,
		n_threshold = 0,
		cor_method = 's',
		cor_pairing = 'pair',
		Forced_zeros = NULL,
		Stratify = NULL
		)
	{
		if(n_b < 1) {
			B_method <- 'Proportion'
			pr_b <- n_b
			n_b <- round(dim(Data)[1] * pr_b, 0)
			if(!n_b) stop('Subsample proportion of data size too small! Must exceed 1/n')
		} else {
			if(n_b%%1!=0) stop('If n_b > 1, n_b must be a whole number.')
			pr_b <- n_b/dim(Data)[1]
			B_method <- 'Static Count'
		}
		
		if(n_b > dim(Data)[1]) stop('Subsample size exceeds the number of records!')
		E_aggr <- list(NA)
		if(is.null(V_set)) V_set <- names(Data)
		k <- length(V_set)
		M <- data.frame(b = numeric(0), L1 = numeric(0), t = numeric(0), m = numeric(0))

		B_sets <- sapply(rep(n_b, B), sample, x = 1:dim(Data)[1])
		Data_b <- array(0, dim = c(n_b, k, B),  dimnames = list(1:n_b, V_set, 1:B))
		W_aggr <- array(0, dim = c(k, k, B, length(L1_set)), dimnames = list(V_set, V_set, 1:B, as.character(L1_set)))

		iter <- 1
		for(b in 1:B) {
			t_0 <- proc.time()[3]
			Data_b[,,b] <- as.matrix(Data[B_sets[,b],V_set])
			for(L1 in L1_set){
				R_i <- Rnet(
					as.data.frame(Data_b[,,b]), 
					L1, 
					cor_method = cor_method, 
					cor_pairing = cor_pairing,
					n_threshold = n_threshold,  
					Forced_zeros = Forced_zeros,
					Stratify = Stratify
					)
				if(ecount(R_i@R) > 0 ) {
					E_i <- as_edgelist(R_i@R)
					E_aggr[[iter]] <- data.frame(
						b = b, 
						L1 = L1, 
						V1 = E_i[,1], 
						V2 = E_i[,2], 
						E = paste(E_i[,1], E_i[,2], sep = '--'),
						omega = E(R_i@R)$omega
						)
					M[iter,] <- c(b, L1, proc.time()[3] - t_0, ecount(R_i@R))
					W_aggr[,,b, as.character(L1)] <- R_i@Omega
					iter <- iter + 1
				} else break()
			}
		}
		E_aggr <- as.data.frame(rbindlist(E_aggr))
		Comp_Stab <- aggregate(Count ~ L1 + E, data = cbind(E_aggr[c('L1', 'E')], Count = 1), FUN = sum)
		Comp_Stab$Pr <- Comp_Stab$Count/B
		Comp_Stab <- within(Comp_Stab, Eta <- 2 * Pr * (1-Pr)/choose(k, 2) )
		StARS <- aggregate(Eta ~ L1, Comp_Stab, FUN = sum)
		StARS.vec <- StARS$Eta
		names(StARS.vec) <- StARS$L1
		print(paste('L1Selection completed', iter, 'loops over', B, 'subsets in', round(sum(M$t), 2), 'seconds'))
				
		return(new('rnet.L1.set',
			Data = Data,
			L1_set = L1_set,
			B = B,
			B_method = B_method,
			B_sets = B_sets,
			Data_b = Data_b,
			pr_b = pr_b,
			n_b = n_b,
			W_aggr = W_aggr,
			E_aggr = E_aggr,
			M = M,
			Edge_stability = Comp_Stab,
			StARS_D = StARS.vec
			))
	})

