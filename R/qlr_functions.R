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
#' @tests
#' example_cor = random_cor_cov(size=5, cors="high", sds=c(1,2,3,4,5))
#' testthat::expect_true(length(names(example_cor))==3)
#' testthat::expect_true(example_cor$sd[1] == 1)
#' example_cor = random_cor_cov(size=5, cors="high")
#' testthat::expect_true(all(example_cor$sd^2 == diag(example_cor$cov.mat)))
random_cor_cov = function(size=6, cors=c("high", "mid", "low"), sds=NULL){

    cor_matrix = fill_matrix(size, cors, check.pd=TRUE)
    
    if (is.null(sds)){
        sds = runif(size, .2, 10)
    }
    
    cov.mat = diag(sds)%*%cor_matrix%*%t(diag(sds))
    
    
    list(cor = cor_matrix, sd= sds, cov.mat = cov.mat)
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
simulate_studies = function(rho_matrix, number_of_studies=10, prob_any_missing = .5, prop_missing = .6, imputations=5){
    # simulate for a bunch of studies
    num_variables = nrow(rho_matrix$cor)
    study_correlations = 1:number_of_studies %>% map(map_studies,
                                                     rho = rho_matrix,
                                                     prob_any_missing=prob_any_missing,
                                                     prop_missing=prop_missing)
    
    # input missing correlations and convert to columns
    study_correlations_NAd = study_correlations %>%
        map_dfc(return_na_correlations, letters[1:num_variables]) %>%
        t %>%
        data.frame() 
    names(study_correlations_NAd) = name_vechs(letters[1:num_variables], collapse = "_") 
    study_correlations_NAd = fill_missing_columns(study_correlations_NAd)
    # impute the missing values
    #fill_missing_columns(study_correlations_NAd)
    imputed_matrix = mice::mice(data.matrix(study_correlations_NAd), m=imputations, remove.collinear=FALSE)
    imputed_summaries = summarize_imputations(imputations, imputed_matrix)
    return(imputed_summaries)
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

summarize_imputations = function(imps, imputed_matrix, vechs=T) {
    imputed_summary = 1:imps %>% 
        map_df(summarize_single_imputation, imputed_matrix) %>% 
        group_by(set) %>% 
        summarize(mean_mean = mean(M),
                  vb = sd(M),
                  vw = var(S)) %>% 
        mutate(pooled_variance = vw + vb + vb/imps) %>% 
        select(set, mean_mean, pooled_variance) %>% 
        set_names(c("Parameter", "Pooled_Estimate", "Pooled_Variance"))
    return(imputed_summary)
    
}

#' Return the limits of a correlation
#'
#' @param cors The size of the correlations between variables. Current options
#' include "high", "mid", and "low". Future iterations may include "toeplitz" and maybe
#' other options
#'
#' @return A vector containing the ranges of the correlations
#' @tests 
#' testthat::expect_equal(return_cor_limits("high")[1], 0.5)
#' testthat::expect_equal(return_cor_limits("mid")[1], 0.2)
#' testthat::expect_equal(return_cor_limits("low")[1], 0.05)
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
#' @tests 
#' example_cor = fill_matrix(size=5, cors="high")
#' testthat::expect_true(all(example_cor>.5))
#' testthat::expect_true(det(example_cor)>0)

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



#' @tests
#' testthat::expect_true(isSymmetric(make_matrix_symmetric(diag(1, nrow=3), runif(3))))
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

#' @tests
#' cor_mat = fill_matrix(10,cors = "mid")
#' data = MASS::mvrnorm(100, mu=rep(0, times=10), Sigma = cor_mat)
#' testthat::expect_true(all(head(data)==head(remove_some_variables(data, prob_any_missing = 0))))
#' testthat::expect_true(ncol(remove_some_variables(data, 1, .2))==8)
remove_some_variables = function(data, prob_any_missing = 1, proportion_missing=0) {
    
    if (runif(1)>prob_any_missing) return(data)
    
    total_rows = ncol(data)
    num_missing = round(total_rows*proportion_missing)
    
    ## sample which ones are removed
    rows_to_remove = sample(1:total_rows, size=num_missing)
    
    return(data[, -rows_to_remove])    
}

#' @tests
#' # n param either TRUE or a number. If TRUE, it will randomly sample a sample size. 
#' # If a number, it will sample with that n
#' testthat::expect_true(isSymmetric(map_studies(1, random_cor_cov(cors="low"), n=100)))
map_studies = function(number, rho, prob_any_missing=.5, prop_missing=.2, n=T) {

    if (n == TRUE) n = sample(c(30, 100, 200, 500, 1000), size=1)
    return(simulate_a_study(rho, n, prob_any_missing, prop_missing))
    
}

#' @importFrom OpenMx vechs
#' @tests
#' # a function that takes a correlation, takes ALL variables, and puts NA where variables are missing
#' rho = random_cor_cov(size=5, cors="high")
#' cor_mat = simulate_a_study(rho, n=100, prob_any_missing=1, proportion_missing = .6)
#' variable_names = letters[1:5]
#' length_varnames = length(variable_names)
#' testthat::expect_true(ncol(return_na_correlations(cor_mat, variable_names, vechs=FALSE)) == length_varnames)
#' testthat::expect_true(length(return_na_correlations(cor_mat, variable_names, vechs=TRUE)) == (length_varnames*(length_varnames-1))/2)
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

#' @tests
#' testthat::expect_true(name_vechs(letters[1:3])[3]=="b:c")
name_vechs = function(variable_names, collapse=":"){
    combn(variable_names, 
          m=2, 
          FUN = function(x) paste0(x[1], collapse, x[2]))
}


# correlations_dframe = study_correlations_NAd
# rho = random_cor_cov(size=5, cors="high")
# cor_mat = simulate_a_study(rho, n=100, prob_any_missing=1, proportion_missing = .8)
fill_missing_columns = function(correlations_dframe){
    #browser()
    all_missing = apply(correlations_dframe, 2, function(x) all(is.na(x)))
    if (!all(all_missing)) return(correlations_dframe)

    random_correlations = runif(nrow(correlations_dframe)*length(sum(all_missing)),
                                -1, 1)
    correlations_dframe[,all_missing] = random_correlations
}
 