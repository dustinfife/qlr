## read in the data
d = read.csv("inst/data/schizo_qlr_aug2020.csv", sep=";")
head(d)
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
