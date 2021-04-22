require(tidyverse)

## read in the data
d = read.csv("inst/data/schizo_qlr_aug2020.csv", sep=";", stringsAsFactors = FALSE) %>% 
  filter(estimate_type == "Correlation" & variable_1  != "" & 
           variable_1 != "O-LIFE Total" & variable_2 != "O-LIFE Total") %>% 
  dplyr::select(value, study, article_id, variable_1, variable_2) %>% 
  filter(abs(value)<1)
d$variable_1[d$variable_1=="SPQ-B Positive"] = "SPQ Positive"
d$variable_1[d$variable_1=="SPQ-B Negative"] = "SPQ Negative"
d$variable_1[d$variable_1=="SPQ-B Disorganized"] = "SPQ Disorganized"

d$variable_2[d$variable_2=="SPQ-B Positive"] = "SPQ Positive"
d$variable_2[d$variable_2=="SPQ-B Negative"] = "SPQ Negative"
d$variable_2[d$variable_2=="SPQ-B Disorganized"] = "SPQ Disorganized"

## read in table with the different names and remap them
names = read.csv("inst/data/names_db.csv")

# create function that maps names to qualtrics names
convert_names = function(x) {
  return_value = as.character(names$model_names[which(names$db_names == x)][1])
  if (return_value == "") {
    x = gsub(" ", "_", x) 
    return(x)
  }  
  return(return_value)
}
d = d %>% mutate(
  variable_1 = variable_1 %>% map_chr(convert_names),
  variable_2 = variable_2 %>% map_chr(convert_names)
  )


## convert to list format
listed_correlations = matrix_lists(d, "article_id", v1 = "variable_1", v2="variable_2")

## convert to column correlation format
column_correlations = correlations_to_columns(listed_correlations, fill_missing = F)
# 40 variables


# rename columns (otherwise, mice has issues)
old_names = names(column_correlations)
new_names = paste0("a", 1000:(1000+ncol(column_correlations)-1))
names(column_correlations) = new_names


### set that loop
results_symptoms = as.list(1:1000)
results_measure = as.list(1:1000)

  for (i in 1:length(results_symptoms)){
  
  # find columns that can be imputed with MICE (at least two estimates)
  imputable_columns = apply(column_correlations, 2, function(x) length(which(!(is.na(x))))>1)
  imputed = impute_correlations(column_correlations[,imputable_columns], imputations = 1)
  
  # rename columns
  names(column_correlations) = old_names
  names(imputed)[-1] = old_names[imputable_columns]
  
  # fill in missing correlations and compute means
  imputed_matrix = correlations_to_columns(listed_correlations, fill_missing = T)
  imputed_matrix[,imputable_columns] = imputed[,-1]
  study_j_matrix = as.matrix((unvechs(imputed_matrix[1,])))
  
  # prepare covariance matrix for lavaan
  require(lavaan)
  study_j_matrix = data.matrix(study_j_matrix)
  variable_names = listed_correlations %>% map(colnames) %>% unlist %>% unique
  colnames(study_j_matrix) = variable_names
  rownames(study_j_matrix) = variable_names
  keepers = variable_names %in% names$model_names
  study_j_matrix = study_j_matrix[keepers,keepers]  
  study_j_matrix = data.matrix(Matrix::nearPD(study_j_matrix, corr=TRUE)$mat)
  det(study_j_matrix)
  # now, fit the SEM for each individual study, multiple times
  measure.model ='
  schizotypy =~ spqbr_pos + spqbr_neg + spqbr_dis
  olife =~olife_ue + olife_ia + olife_in + olife_cd
  mssb =~ mssb_pos + mssb_neg + mssb_dis
  chap =~ wss_sa + wss_phys + wss_mi + wss_per
  '
  symptom.model ='
  positive =~ spqbr_pos + mssb_pos + olife_ue + wss_mi + wss_per 
  negative =~ wss_sa + wss_phys + olife_ia + spqbr_neg + mssb_neg
  disorganized =~ spqbr_dis + mssb_dis + olife_cd 
  olife_in ~~ positive
  olife_in ~~ negative
  olife_in ~~ disorganized
  '
  fit_measure = sem(measure.model, sample.cov=study_j_matrix, sample.nobs = 200)
  results_measure[[i]] = fit_measure@Model@GLIST
  
  fit_symptom = sem(symptom.model, sample.cov=study_j_matrix, sample.nobs = 200)
  results_symptoms[[i]] = fit_symptom@Model@GLIST
  
  print(paste0("Iteration ", i, ", of ", length(results_measure), "\n"))

}

tempfunc = function(x, type) {
  m = rowSums(x$lambda)
  names(m) = lavNames(type)
  m
}
trimme = function(x) {
  x[x<(-1)] = NA
  x[x>1] = NA
  x
}

df_measure2 = results_measure[[1:456]]
df_measure = results_measure %>% map_dfr(tempfunc, type=fit_measure) %>% 
  mutate_all(trimme)
head(df_measure)
flexplot(spqbr_dis~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(olife_ia~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(olife_cd~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(mssb_pos~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(mssb_neg~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(mssb_dis~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(wss_sa~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(wss_phys~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(wss_mi~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))
flexplot(wss_per~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))


df_symptoms = results_symptoms %>% map_dfr(tempfunc, type=fit_symptom)
require(flexplot)
flexplot(spqbr_neg~1, data=df_measure) +
  coord_cartesian(xlim=c(-1, 1))


results_symptoms %>% map_dfr(tempfunc, type=fit_symptom)

results_measure[[1]]









# keep looping through until a PD matrix is there that doesn't change the MI values
j = 0
study_j_matrix=22
while (j != study_j_matrix) {
  j = study_j_matrix
  study_j_matrix = as.matrix(Matrix::nearPD(unvechs(imputed_matrix[1,]), corr=TRUE )$mat)
  study_j_matrix[lower.tri(study_j_matrix)[imputable_columns]] 
}


# fit a random effects model
require(lme4)
imputed_matrix$Study = 1:nrow(imputed_matrix)
lmer(`MSS Negative_MSS Positive`~1 + (1|Study), data=imputed_matrix)
imputed_matrix[1:5,1:3]

mean_correlations = colMeans(correlations_to_columns(listed_correlations, fill_missing = T))
mean_correlations[which(imputable_columns)] = imputed$Pooled_Estimate



# put into a fixed effect matrix
namevars = unique(c(as.character(d$variable_1), as.character(d$variable_2)))
mega_matrix = diag(1, nrow=length(namevars), ncol=length(namevars))
mega_matrix[lower.tri(mega_matrix)]= mean_correlations
mega_matrix[upper.tri(mega_matrix)]= mean_correlations
imputed_matrix$Study = NULL

require(metaSEM)
?tssem2













OpenMx::vech2full(mean_correlations)[1:5, 1:5]
hist(unlist(mean_correlations))

typeof((mean_correlations))

d = read.csv("inst/data/schizo_qlr_aug2020.csv", sep=";") %>% 
  filter(estimate_type == "Correlation") %>% 
  filter(abs(value)<1) %>% 
  unite(combined, c(variable_1,variable_2), sep="---") %>% 
  group_by(combined) %>% 
  summarize(n=length(combined), median = median(value), mean = mean(value), sd=sd(value)) %>% 
  filter(n>1) %>% 
  arrange(-n)

head(d, n=20)


require(tidyverse)
## subset to just correlations
d2 = d %>% filter(estimate_type == "Correlation") %>% 
  select(value, study, article_id, variable_1, variable_2) %>% 
  filter(abs(value)<1)
unique(d2$article_id) %>% 
  map()
  
  map(colnames)
  unite(combined, c(variable_1,variable_2), sep="---") %>% 
  group_by(combined) 


%>% 
  summarize(n=length(combined), median = median(value), mean = mean(value), sd=sd(value)) %>% 
  filter(n>1) %>% 
  arrange(-n)
write.csv(cors, file="data/schizo_qlr_summary.csv")


rho_matrix = random_cor_cov(size = 5, cors="mid")
simulate_studies(rho_matrix, 10, prob_any_missing = .1, prop_missing = .4)
