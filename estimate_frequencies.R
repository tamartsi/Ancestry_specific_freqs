# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
estimate_frequencies <- function(allele_counts, prop_mat, confidence = 0.95, 
                                 low_freq_bound = 0.001, high_freq_bound = 0.999,
                                 use_smoothing_data = FALSE){
  stopifnot(length(allele_counts) == nrow(prop_mat))
  
  prop_mat <- as.matrix(prop_mat)
 
  # turn each allele dosage into a two rows, each reporting one allele
  # with the same proportion ancestries.
  decomposed_alleles <- sapply( allele_counts,decompose_two_alleles_one_person)
  decomposed_alleles <- as.numeric(matrix(decomposed_alleles)	)
  prop_mat_double <- prop_mat[rep(1:nrow(prop_mat), each = 2),]
  
  # add made-up data to avoid boundaries of the frequency parameter space   
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
    prop_mat <- rbind(prop_mat_double, smoothing_data$simulated_prop_mat)
    allele_counts <- c(decomposed_alleles, smoothing_data$simulated_allele_count)
  }
  
  
  ## compute the negative log likelihood function
  nll <- function(freqs){
    allele_probs <- as.numeric(prop_mat_double %*% freqs)
    nll_by_obs <- decomposed_alleles*log(allele_probs) + 
      (1-decomposed_alleles)*log(1-allele_probs)
    return(-sum(nll_by_obs))		
  }
  fit <- optim(par = rep(1/ncol(prop_mat), ncol(prop_mat)), fn = nll, hessian = TRUE, 
               lower = rep(low_freq_bound ,ncol(prop_mat)), upper = rep(high_freq_bound, ncol(prop_mat)), 
               method = "L-BFGS-B")
  
  estimated_freqs <- fit$par
  hessian <- fit$hessian
  ses <- sqrt(diag(solve(hessian)))
  res <- data.frame(ancestry = colnames(prop_mat), 
                    estimated_freq = estimated_freqs, 
                    low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                    high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
  
  return(res)
}
