require(tidyverse)
require(OpenMx)
require(flexplot)

num_variables = 6
size_correlation = sample(c("low", "mid", "high"), size=1)
number_of_studies = 10


display_all = function(i, rho_matrix, mc_final) {
  pop = vechs(rho_matrix$cor)
  val = unique(mc_final$Parameter)[i]
  new_d = mc_final[mc_final$Parameter==val,]
  return(flexplot(Pooled_Estimate~1, data=new_d) + geom_vline(xintercept = pop[i], col="red") + coord_cartesian(xlim=c(0, 1)))
}


# create population correlations
set.seed(23443)
rho_matrix = random_cor_cov(size = num_variables, cors=size_correlation)
f = function(i, rho_matrix) {print(i); simulate_studies(rho_matrix, 10, .8, .7)}
mc_results = 1:200 %>% map(safely(f), rho_matrix)
mc_final = mc_results %>% map(pluck("result")) %>% bind_rows(.id = "iteration") 
#save(mc_final, file="inst/data/mc_v2.Rdata")
mc_final %>% dplyr::filter(Parameter == "a_b") %>% flexplot(~Pooled_Estimate~1, data=.) + geom_vline(xintercept = rho_matrix$cor[1,2])
a = 1:length(unique(mc_final$Parameter)) %>%  map(display_all, rho_matrix, mc_final)
cowplot::plot_grid(plotlist=a)

