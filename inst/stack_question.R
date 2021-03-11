# this function takescolumn data matrix and turns into list ofmatrices
columns_to_lists = function(group_correlations) {
  
  # find names of all variables
  unique_vars = unique(c(group_correlations$v1, group_correlations$v2))
  all_names = combn(unique_vars, 2, simplify=F)
  
  # create empty  matrix of data
  ordered_values = map_dbl(all_names, function(x) assign_correlations(x, group_correlations))
  cor_mat = diag(nrow=length(unique_vars), ncol=length(unique_vars))
  cor_mat[lower.tri(cor_mat)] = ordered_values
  cor_mat[upper.tri(cor_mat)] = ordered_values
  cor_mat = data.frame(cor_mat)
  names(cor_mat) = rownames(cor_mat)= unique_vars
  return(cor_mat)
}

assign_correlations = function(x, group_correlations) {
  #browser()
  id = which(group_correlations$v1 == x[1] & group_correlations$v2 == x[2])
  if (length(id)==0) return(NA)
  return(group_correlations$value[id][1])
}

d = data.frame(id=sort(rep(1:5, times=5)), 
               value = runif(25), 
               v1=sample(letters[1:4], 25, replace=T), 
               v2=sample(letters[1:4], 25, replace=T))
d %>% group_by(id) %>% group_split() %>% 
  map(columns_to_lists)



group_correlations = data.frame(id = c(4,4,4), value = c(.1, .2, .3), v1 = c("a", "a", "b"), v2 = c("b", "c", "c"))
group_correlations %>% arrange(v1, v2)
columns_to_lists(group_correlations)



assign_correlations("a", "b", group_correlations)


vech2full(group_correlations$value)
vech2full(1:10)

matrix(1:16, 4, 4)
vech(matrix(1:16, 4, 4))
vech2full(vech(matrix(1:16, 4, 4)))


x <- list(1, 1, 1)
y <- list(10, 20, 30)
z <- list(100, 200, 300)

map2(x, y, ~ .x + .y)
x = letters[1:3]
y = letters[4:6]
map2(x,y, function(a,b) paste0(a, "-", b))



?combn
