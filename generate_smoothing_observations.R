# Generate artificial observations with alleles from each provided ancestries.
# These observations are later added to the real data to protect from reaching boundary condition
# and failure of the estimation function. 

.generate_smoothing_observations <- function(ancestry_names){
  
  n_ancestry <- length(ancestry_names)
  
  # construct a matrix with two rows per ancestry, corresponding
  # to two observations with 100% global proportions of that ancestry
  prop_mat <- matrix(0, nrow =n_ancestry*2, ncol = n_ancestry)
  colnames(prop_mat) <- ancestry_names
  
  for (i in 1:n_ancestry){
    row_inds <- (i-1)*2+c(1,2)
    prop_mat[row_inds,i] <- 1
  }
   
  # for each ancestry, we set one observation with one variant allele
  # and a second observations with non variant alleles. 
  allele_counts <- rep(c(1,0), n_ancestry)
  
  return(list(simulated_prop_mat = prop_mat, simulated_allele_count = allele_counts))
   
}