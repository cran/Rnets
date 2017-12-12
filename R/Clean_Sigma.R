#' Clean_Sigma function
#'
#' Primarily a convenience function, Clean_Sigma estimates a correlation (Sigma) matrix that may contain substantial amounts of missing data and returns it in an easily usable form for the glasso function.
#' @param x Rhe data used to estimate the correlation matrix. Pairwise correlations between all columns are estimated.
#' @param cor_method Method used to estimate correlations. Defaults to Spearman's method. See help file for 'cor' function.
#' @param cor_pairing Method used to determine how NA entries are handled. See help file for 'cor' function.
#' @importFrom stats cor
#' @export

Clean_Sigma <- function(x,
				cor_method = 's',
				cor_pairing = 'pair')
{
	CorrMat <- 	cor(x, 
				method = cor_method, 
				use = cor_pairing
				)
	CorrMat[is.na(CorrMat)] <- 0
	diag(CorrMat) <- 1
	return(CorrMat)
}