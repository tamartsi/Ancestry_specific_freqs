# Generate artificial observations with alleles from each provided ancestries.
# These observations are later added to the real data to protect from reaching boundary condition
# and failure of the estimation function. 

.generate_smoothing_observations_lafa <- function(ancestry_names){
  
  ancestry1 <- rep(ancestry_names, each = 2)
  allele_counts <- rep(c(1,0), length(ancestry_names))
  
  return(list(simulated_ancestry1 = ancestry1, simulated_allele_count = allele_counts))
   
}
