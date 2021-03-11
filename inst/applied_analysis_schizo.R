require(tidyverse)

## read in the data
d = read.csv("inst/data/schizo_qlr_aug2020.csv", sep=";") %>% 
  filter(estimate_type == "Correlation" & variable_1  != "") %>% 
  select(value, study, article_id, variable_1, variable_2) %>% 
  filter(abs(value)<1)

## convert to list format
listed_correlations = matrix_lists(d, "article_id", v1 = "variable_1", v2="variable_2")

## convert to column correlation format
column_correlations = correlations_to_columns(listed_correlations, fill_missing = F)

# rename columns (otherwise, mice has issues)
old_names = names(column_correlations)
new_names = paste0("a", 1000:(1000+ncol(column_correlations)-1))
names(column_correlations) = new_names

# find columns where all are missing, then remove from imputation
all_missing = apply(column_correlations, 2, function(x) all(is.na(x)))
column_correlations[,!all_missing][,1:5]
imputed = impute_correlations(column_correlations[,!all_missing])
## imputations predict NA for situations where there'sonly one estimate
## and thecolumns aren'tlining up

names(column_correlations) = old_names
imputed$Parameter = old_names[!all_missing]

d[(d$variable_1=="MSS Negative" & d$variable_2=="MSS Disorganized") | 
    (d$variable_2=="MSS Negative" & d$variable_1=="MSS Disorganized"),]


# compute the means
mean_correlations = colMeans(correlations_to_columns(listed_correlations, fill_missing = T))
mean_correlations[which(!all_missing)] = imputed$Pooled_Estimate

column_correlations[1:5,which(!all_missing)[1:5]]
mean_correlations[which(!all_missing)[1:5]]

# put into a matrix
namevars = unique(as.character(d$variable_1), as.character(d$variable_2))
mega_matrix = diag(1, nrow=length(namevars), ncol=length(namevars))
mega_matrix[lower.tri(mega_matrix)]= mean_correlations
dim(mega_matrix)
length(mean_correlations)
length(lower.tri(mega_matrix))
(136*135)/2
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
