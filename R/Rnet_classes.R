#' An S4 class for accepting input data common for generating all rnet objects. These objects are used to handle rnet objects and need not be called by the user.
#' 
#' @slot RawData A dataframe containing the original dataset.
#' @slot cor_method The type of correlation matrix to use. Must be a partial match to one of the following strings: 'pearson', 'spearman', 'kendall'.
#' @slot cor_pairing Method for handling how missing data is handled in pairs. See 'pair' argument in function 'cor' for more information.
#' @slot n_threshold The minimum number of valid pair-wise observations that must exist for an edge to be estimated. Vertex pairs with fewer valid pair-wise observations are assumed to be conditiontally independent.
#' @slot L1_orig The declared L1 penalty to be used when estimating rnet topology
#' @slot v_set_orig The declared set of k variables to be included in the rnet as vertices
#' @slot Forced_zeros A matrix with 2 columns containing pairs of vertices to force to be conditionally independent in the rnet
#' @slot Layout_master A k x 2 matrix x & y coordinates of each vertex in the graph.
#' @rdname rnetInput

rnetInput <- setClass(Class = "rnetInput",
	slots = c(
		RawData = 'data.frame',
		cor_method = 'character',
		cor_pairing = 'character',
		n_threshold = 'numeric',
		L1_orig = 'numeric',
		V_set_orig = 'character',
		Forced_zeros = 'matrix',
		Layout_master = 'ANY'
		)
	)

#' An S4 class containg the information of a basic rnet. These objects are used to handle rnet objects and need not be called by the user.
#' 
#' Inherits from 'rnetInput'
#' @slot Data A dataframe containing the dataset.
#' @slot L1 The L1 penalty
#' @slot V_set The vertex set
#' @slot n A square matrix of n_ij, the number of valid pairs used to estimate the edge between i and j.
#' @slot V_ommitted The names of vertices removed from the rnet due to low sample size or forced to 0.
#' @slot V_metadata A Vector containing the attribute names assigned to the rnet's vertices.
#' @slot E_metadata A vector containing the attribute names assinged to the rnet's edges.
#' @slot Zeros A list of edges forced to zero due to small sample size or declared by user.
#' @slot Sigma The empirical correlation matrix
#' @slot Theta The penalized precision matrix (Sigma^-1)
#' @slot Omega The penalized partial correlation matrix. Also the weighted adajacency matrix used for the network.
#' @slot A The adjacency matrix (omega_ij != 0)
#' @slot R An igraph object derived from Omega
#' @slot Layout A k' x 2 matrix containing the x & y coords for the vertices for plotting. k' indicates that some of the k vertices declared may have been removed.
#' @rdname rnetBasic
#' 
rnetBasic <- setClass(Class = "rnetBasic",
	slots = list(
		Data = 'data.frame',
		L1 = 'numeric',
		V_set = 'character',
		n = 'matrix',
		V_omitted = 'character',
		V_metadata = 'character',
		E_metadata = 'character',
		Zeros = 'list',
		Sigma = 'matrix',
		Theta = 'matrix',
		Omega = 'matrix',
		A = 'matrix',
		loglik = 'numeric',
		R = 'ANY',
		Layout = 'matrix'
		),
	contains = 'rnetInput'
	)

#' An S4 class containing the information of an rnet representing one stratum of data. These objects are used to handle rnet objects and need not be called by the user.
#'
#' Inheirits from 'rnetBasic'
#' @slot Strata_def The expression used to define the stratum of data represented by the rnet, e.g., if Statra_def = as.expression(Year == 2008), the rnet is estimated from the subset of data collected in the Year 2008.
#' @rdname rnetStrata
rnetStrata <- setClass(Class = "rnetStrata",
	slots = list(Strata_def = 'expression'),
	contains = 'rnetBasic'
	)
#' An S4 class containing multiple Rnets, each the information of an rnet representing one stratum of data. These objects are used to handle rnet objects and need not be called by the user.
#'
#' Inheirits from 'rnetInput'
#' @slot Stratify_by The name of the variable used to the original data, e.g. if Stratify_by = 'Year', the object contains 1 rnetStrata object for each unique value of the variable "year" in the dataset.
#' @slot E_matrix A matrix containing every edge found in the entire set of rnets, and the stratum in which it was found.
#' @slot R_strata A list of rnetStrata objects for each strata of the declared Stratify_by variable.
#' @rdname rnetMultiStrata

rnetMultiStrata <- setClass(Class = "rnetMultiStrata",
	slots = list(
		Stratify_by = 'character',
		E_matrix = 'matrix',
		R_Strata = 'list'
		),
	contains = 'rnetInput'
	)

#' An S4 class for tracking multiple rnets created from one data set over multiple L1 penalties for comparison and selection.
#' 
#' See "L1_selection" documentation for more information.
#' @slot Data A dataframe containing the dataset
#' @slot L1_set a numeric vector containing the candidate L1 penalties
#' @slot B The number of subsamples to draw from the data to evaluate topologic stability 
#' @slot B_method Assigned either "proportionate" or "Total number" depending on how subsample size is determined.
#' @slot B_sets A matrix (n_B x B) containing the rownumbers of the subsamples
#' @slot Data_b An arrary (n_b x k x B) containing the data from the B_sets matrix
#' @slot pr_b The size of the subsample B as a proportion of the complete dataset
#' @slot n_b The size of the subsample
#' @slot W_aggr An array (k x k x B x L1) containg all the weighted adjacency matrices generated by the all of subsamples over the L1 penalties
#' @slot A_aggr An array (k x k x B x L1) containg all the adjacency matrices generated by the all of subsamples over the L1 penalties
#' @slot M A dataframe with with graphical density data over the set of generated networks.
#' @slot Edge_stability A dataframe showing edge stability over the set of generated networks.
#' @slot StARS_D A vector of D_b values used for L1 selection.
#' @rdname rnetInput
#' 
rnet.L1.set <- setClass(Class = 'rnet.L1.set',
	slots = list(
		Data = 'data.frame',
		L1_set = 'numeric',
		B = 'numeric',
		B_method = 'character',
		B_sets = 'matrix',
		Data_b = 'array',
		pr_b = 'numeric',
		n_b = 'numeric',
		W_aggr = 'array',
		E_aggr = 'data.frame',
		M = 'data.frame',
		Edge_stability = 'data.frame',
		StARS_D = 'numeric'
		)
	)