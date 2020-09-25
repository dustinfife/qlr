# source("R/qlr_functions.R")
# require(tidyverse)
# require(OpenMx)
# 
# num_variables = 6
# size_correlation = sample(c("low", "mid", "high"), size=1)
# number_of_studies = 10
# 
# # create population correlations
# rho_matrix = random_cor_cov(size = num_variables, cors=size_correlation)
# 
# # simulate for a bunch of studies
# study_correlations = 1:number_of_studies %>% map(map_studies, 
#                             rho = rho_matrix, 
#                             prob_any_missing=.8, 
#                             prop_missing=.6)
# 
# # input missing correlations and convert to columns
# study_correlations_NAd = study_correlations %>% 
#     map_dfc(return_na_correlations, 
#             letters[1:num_variables]) %>% 
#     t %>% 
#     data.frame() 
# names(study_correlations_NAd) = name_vechs(letters[1:6], collapse = "_")
# 
# 
# 
# # impute the missing values
# require(mice)
# typeof(study_correlations_NAd)
# imputed_matrix = mice(data.matrix(study_correlations_NAd), m=5, method="rf")
# complete(imputed_matrix, 1)
# study_correlations_NAd
# 
# summarize_imputation = function(i)
# means = 1:5 %>% map()
# 
# #### loop through and average the correlations
# mean.cors = data.frame(matrix(nrow=imps, ncol=unique.cors))
# sd.cors = data.frame(matrix(nrow=imps, ncol=unique.cors))
# i = 1
# for (i in 1:nrow(mean.cors)){
#     newd = complete(id, i)
#     mean.cors[i,] = colMeans(newd[,-1])
#     sd.cors[i,] = apply(newd[,-1], 2, sd)
# }
# colMeans(mean.cors)[c(1,2,p)]
# apply(mean.cors,2,sd)[c(1,2,10)] + colMeans(sd.cors)[c(1,2,10)]
