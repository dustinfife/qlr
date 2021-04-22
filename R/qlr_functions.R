#' Generate a random correlation (or covariance) matrix
#'
#' @param size The number of variables. Defaults to 6.  
#' @param cors The size of the correlations between variables. Current options
#' include "high", "mid", and "low". Future iterations may include "toeplitz" and maybe
#' other options
#' @param sds A vector containing the standard deviation of each of the variables. 
#' Defaults to NULL, in which case, it randomly generates standard deviations that range
#' from .2 to 10
#'
#' @return A list containing
#' \itemize{
#'   \item cor - the randomly generated correlation matrix
#'   \item sd - a vector of standard deviations
#'   \item cov.mat - the variance/covariance matrix
#' }
#' @export
#'
#' @examples
#' random_cor_cov(size=4, cors="low")
#' random_cor_cov(size=3, cors="mid", sds = c(5,5,5))
random_cor_cov = function(size=6, cors=c("high", "mid", "low"), sds=NULL){

    cor_matrix = fill_matrix(size, cors, check.pd=TRUE)
    
    if (is.null(sds)){
        sds = runif(size, .2, 10)
    }
    
    cov.mat = diag(sds)%*%cor_matrix%*%t(diag(sds))
    
    
    list(cor = cor_matrix, sd= sds, cov.mat = cov.mat)
}

# this function takescolumn data matrix and turns into list ofmatrices
#' @importFrom dplyr group_by group_split
#' @importFrom rlang sym `!!`
#' @importFrom purrr map
matrix_lists = function(d, id="id", v1="v1", v2="v2", value="value") {
    
    # make sure all variables are in the dataset
    user_supplied_names = c(id, v1, v2, value)
    if (!(all(user_supplied_names %in% names(d)))){
        missingvar = which(!(user_supplied_names %in% names(d)))
        stop(paste0(user_supplied_names[missingvar], 
                    " is not a variable in your dataset"))
    }
    
    # loop through and convert to list
    id_name = rlang::sym(id)
    d = d %>% dplyr::group_by(!!(id_name)) %>% dplyr::group_split() %>% 
        purrr::map(columns_to_lists, v1=v1, v2=v2, value=value)
    return(d)
}

    #' @importFrom purrr map_dbl
    columns_to_lists = function(group_correlations, v1="v1", v2="v2", value="value") {
        
        group_correlations[[v1]] = as.character(group_correlations[[v1]])
        group_correlations[[v2]] = as.character(group_correlations[[v2]])
        # find names of all variables
        unique_vars = unique(c(group_correlations[[v1]], group_correlations[[v2]]))
        all_names = combn(unique_vars, 2, simplify=F)
        
        # create empty  matrix of data
        ordered_values = purrr::map_dbl(all_names, 
                function(x) assign_correlations(x, group_correlations, v1=v1, v2=v2, value=value))
        cor_mat = diag(nrow=length(unique_vars), ncol=length(unique_vars))
        cor_mat[lower.tri(cor_mat)] = ordered_values
        cor_mat[upper.tri(cor_mat)] = ordered_values
        cor_mat = data.frame(cor_mat)
        names(cor_mat) = rownames(cor_mat)= unique_vars
        return(cor_mat)
    }
    #helper function for columns_to_lists
    # given a pair of variable names, this will return the corresponding value
    assign_correlations = function(x, group_correlations, v1="v1", v2="v2", value="value") {
      #
      id = which(group_correlations[[v1]] == x[1] & group_correlations[[v2]] == x[2])
      if (length(id)==0) return(NA)
      return(group_correlations[[value]][id][1])
    }

#' Simulate a chosen number of studies based on a fixed correlation matrix
#'
#' @param rho_matrix The correlation matrix (obtained from \link{random_cor_cov}). 
#' @param number_of_studies The total number of studies to simulate
#' @param prob_any_missing The probability that any given study will have at least one missing correlation
#' @param prop_missing The proportion of variables from each study that will have missing data
#'
#' @return
#' @export
#' @importFrom mice mice
simulate_studies = function(rho_matrix, number_of_studies=10, prob_any_missing = .5, prop_missing = .6, imputations=5, return.list=FALSE){

    # simulate for a bunch of studies
    num_variables = nrow(rho_matrix$cor)
    study_correlations = 1:number_of_studies %>% map(map_studies,
                                                     rho = rho_matrix,
                                                     prob_any_missing=prob_any_missing,
                                                     prop_missing=prop_missing)
    
    if (return.list) return(study_correlations)
    # input missing correlations and convert to columns
    study_correlations_NAd = correlations_to_columns(study_correlations, T)
    
    # impute the missing values
    return(impute_correlations(study_correlations_NAd, imputations))
}

  correlations_to_columns = function(study_correlations, fill_missing = TRUE) {
      
      # extract the variable names
      unique_varnames = study_correlations %>% map(colnames) %>% unlist %>% unique
      matrix_correlations = study_correlations %>%
          map_dfc(return_na_correlations, unique_varnames) %>%
          t %>%
          data.frame() 
      
      # figure out unique variable names
      unique_varnames = study_correlations %>% map(colnames) %>% unlist %>% unique
      names(matrix_correlations) = name_vechs(unique_varnames, collapse = "_") 
      if (fill_missing) return(fill_missing_columns(matrix_correlations))
      return(matrix_correlations)
  }
    
      #' @importFrom OpenMx vechs
      return_na_correlations = function(cor_mat, variable_names, vechs=TRUE) {
  
        ## create NA matrix
        na_cor_mat = data.frame(matrix(NA, nrow=length(variable_names), ncol=length(variable_names)))
        names(na_cor_mat) = row.names(na_cor_mat) = variable_names
        not_missing = row.names(cor_mat)
        
        ## replace NA matrix where we have data
        na_cor_mat[not_missing, not_missing] = cor_mat
        
        if (vechs) return(OpenMx::vechs(na_cor_mat))
        na_cor_mat
      }

  impute_correlations = function(study_correlations, imputations=5, filter_na=TRUE) {
      
      imputed_matrix = mice::mice(data.matrix(study_correlations), 
                                  m=imputations, remove.collinear=FALSE,
                                  method="rf")
      imputed_summaries = summarize_imputations(imputations, imputed_matrix)
      return(imputed_summaries)
  }


unvechs = function(correlations) {
  k = length(correlations)
  num_vars = (1 + sqrt(1 + 8*k))/2
  cor_matrix = diag(1, nrow=num_vars, ncol=num_vars)
  cor_matrix[lower.tri(cor_matrix)] = as.numeric(correlations)
  cor_matrix[upper.tri(cor_matrix)]= as.numeric(correlations)
  cor_matrix
}

summarize_imputations = function(imps, imputed_matrix, vechs=T) {
    imputed_summary = 1:imps %>% 
        # creates a list of data frames
        map(~complete(imputed_matrix, .x)) %>% 
        # gets to data frame level
        map_dfr(~ .x %>% mutate(study = 1:nrow(.x))) %>% 
        group_by(study) %>% 
        summarize_all(mean)
    return(imputed_summary)
}

    summarize_single_imputation = function(i, imputed_matrix) {
      
      complete(imputed_matrix, i) %>% 
        summarise_all(.funs = list(M = ~   mean(x = ., na.rm=T),
                                   S   = ~   var(x = ., na.rm=T))) %>% 
        pivot_longer(everything(),
                     names_to = c("set", ".value"),
                     names_pattern = "(.+)_(.+)"
        )
    }

#' Return the limits of a correlation
#'
#' @param cors The size of the correlations between variables. Current options
#' include "high", "mid", and "low". Future iterations may include "toeplitz" and maybe
#' other options
#'
#' @return A vector containing the ranges of the correlations
return_cor_limits = function(cors=c("high", "mid", "low")) {
    if (cors=="high") return(c(.5, .8))
    if (cors=="mid") return(c(.2, .5))
    if (cors=="low") return(c(.05, .2))
    
    ## return something random if they don't give anything
    runif(2)
    
}



#' Create a correlation matrix, checking for positive definite
#'
#' @param size Number of variables in the variance/covariance matrix
#' @param cors The size of the correlations between variables. Current options
#' include "high", "mid", and "low". Future iterations may include "toeplitz" and maybe
#' other options
#' @param check.pd Should the matrix be checked for a positive definite matrix? Defaults to TRUE. 
#' Probably should be FALSE if the number of variables is high. 
#'
#' @return A randomly-generated correlation matrix
#' @export
#'
#' @examples
#' fill_matrix(size=5, cors="high")
#' fill_matrix(size=25, cors="high", check.pd=FALSE)
fill_matrix = function(size=6, cors=c("high", "mid", "low"), check.pd = TRUE){
    
    cors = match.arg(cors)
    num.cors = (size*(size-1))/2
    cor.mat = diag(1, nrow=size, ncol=size)
    limits = return_cor_limits(cors)
    cor.vals = runif(num.cors, limits[1], limits[2])
    cor_matrix = make_matrix_symmetric(cor.mat, cor.vals)
    
    # return matrix if PD (or if we don't care about PD)
    if (!check.pd | det(cor.mat)>0) return(cor_matrix)
    
    determ = -1
    i = 0
    while (determ < 0 & i<100){
        i = i+1
        
        cor_matrix = make_matrix_symmetric(cor.mat, cor.vals)
        determin = det(cor_matrix)
        if (determin>0) return(cor_matrix)
    }
    
    warning("I couldn't find a positive-definite matrix. I'm returning one that's not positive definite.")
    return(cor_matrix)
}


make_matrix_symmetric = function(cor.mat, cor.vals){
    cor.mat[lower.tri(cor.mat)] = cor.vals
    cor.mat[upper.tri(cor.mat)] = cor.vals
    return(cor.mat)
}

# rho = random_cor_cov(size=5, cors="high")
# study_cor = simulate_a_study(rho=rho, n=30, prob_any_missing=1, proportion_missing=.2)
# testthat::expect_true(all(row.names(study_cor) %in% letters[1:5]))
# testthat::expect_true(nrow(study_cor)==4)        
simulate_a_study = function(rho, n, prob_any_missing, proportion_missing=0, cov=FALSE){
    
    ## simulate data
    data = data.frame(MASS::mvrnorm(n, mu=rep(0, times=nrow(rho$cov.mat)), Sigma=rho$cov.mat))
    names(data) = letters[1:ncol(data)]
    
    ### return original dataset if probabilities work out
    data = remove_some_variables(data, prob_any_missing, proportion_missing=proportion_missing)
    
    ### return correlation
    if (cov) return(cov(data))
    return(cor(data))
}

remove_some_variables = function(data, prob_any_missing = 1, proportion_missing=0) {
    
    if (runif(1)>prob_any_missing) return(data)
    
    total_rows = ncol(data)
    num_missing = round(total_rows*proportion_missing)
    
    ## sample which ones are removed
    rows_to_remove = sample(1:total_rows, size=num_missing)
    
    return(data[, -rows_to_remove])    
}

map_studies = function(number, rho, prob_any_missing=.5, prop_missing=.2, n=T) {
#
    if (n == TRUE) n = sample(c(30, 100, 200, 500, 1000), size=1)
    return(simulate_a_study(rho, n, prob_any_missing, prop_missing))
    
}

name_vechs = function(variable_names, collapse=":"){
    combn(variable_names, 
          m=2, 
          FUN = function(x) paste0(x[1], collapse, x[2]))
}


# correlations_dframe = study_correlations_NAd
# rho = random_cor_cov(size=5, cors="high")
# cor_mat = simulate_a_study(rho, n=100, prob_any_missing=1, proportion_missing = .8)
fill_missing_columns = function(correlations_dframe){
  
    # if there are no missing columns, just return the original dataset
    imputable_columns = apply(correlations_dframe, 2, function(x) length(which(!(is.na(x))))>1)
    if (!any(imputable_columns)) return(correlations_dframe)
    
    # subset the correlations that need to be noninformatively imputed and fill in with random missing values
    new_correlations = as.matrix(correlations_dframe[,!imputable_columns])
    na_cells = is.na(new_correlations)

    # generate with fisher's z
    #random_z = rnorm(length(which(na_cells)), .1)
    #random_r = tanh(random_z)
                 #hist(random_r)               
    random_r = (rbeta(length(which(na_cells)), 4,4)*2-1)
    new_correlations[na_cells] = random_r
    # merge the new imputed data with the old dataset
    correlations_dframe[,!imputable_columns] = new_correlations
    return(correlations_dframe)
}
 