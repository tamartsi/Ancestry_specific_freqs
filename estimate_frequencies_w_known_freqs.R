# Function to estimate ancestry-specific allele frequencies for a given variant,
# when some (one or more) ancestral allele frequencies are known and provided.
# takes a vector of allele_countss allele counts for a given variant for n people 
# and an n x d matrix matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# known_freq need to be a named vector, specifying known frequencies. The names of elements 
# in known_freq are some of the column names of prop_mat. 
estimate_frequencies_w_known_freqs <- function(allele_counts, prop_mat, known_freqs, confidence = 0.95, 
                                            low_freq_bound = 0.001, high_freq_bound = 0.999,
                                            use_smoothing_data = FALSE){
  
  stopifnot(length(allele_counts) == nrow(prop_mat))
  stopifnot(all(is.element(names(known_freqs), colnames(prop_mat))))
  stopifnot(length(known_freqs) > 0 & length(known_freqs) < ncol(prop_mat))

  # turn each allele dosage into a two rows, each reporting one allele
  # with the same proportion ancestries.  
  decomposed_alleles <- sapply( allele_counts,decompose_two_alleles_one_person)
  decomposed_alleles <- as.numeric(matrix(decomposed_alleles)	)
  
  
  # re-order the column so that the ones with known frequencies are at the beginning:
  n_known <- length(known_freqs)
  n_unknown <- ncol(prop_mat) - n_known
  inds_known <- which(is.element(colnames(prop_mat), names(known_freqs)))
  inds_unknown <- setdiff(1:ncol(prop_mat), inds_known)
  prop_mat <- prop_mat[,c(inds_known, inds_unknown)]
  prop_mat_double <- prop_mat[rep(1:nrow(prop_mat), each = 2),]

  # add made-up data to avoid boundaries of the frequency parameter space 
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
    # add only the rows corresponding to the ancestries with unknown frequencies
    inds_keep <- (n_known*2 + 1):((n_known + n_unknown)*2)
    prop_mat <- rbind(prop_mat_double, smoothing_data$simulated_prop_mat[inds_keep,])
    allele_counts <- c(decomposed_alleles, smoothing_data$simulated_allele_count[inds_keep])
  }
  
  
  
  known_sums <- prop_mat_double[,names(known_freqs), drop = FALSE] %*% known_freqs
  
  
  
  ## compute the negative log likelihood function
  nll <- function(unknown_freqs){
    allele_probs <- as.numeric(known_sums + prop_mat_double[,c(n_known+1):ncol(prop_mat_double), drop = FALSE] %*% unknown_freqs)
    nll_by_obs <- decomposed_alleles*log(allele_probs) + 
      (1-decomposed_alleles)*log(1-allele_probs)
    return(-sum(nll_by_obs))		
  }
  if (n_unknown == 1){
    fit <- optim(par = rep(0.5, n_unknown), fn = nll, hessian = TRUE,
                 method = "Brent", lower = low_freq_bound, upper = high_freq_bound)
  } else{
    fit <- optim(par = rep(0.5, n_unknown), fn = nll, hessian = TRUE,
                 lower = rep(low_freq_bound, n_unknown), 
                 upper = rep(high_freq_bound, n_unknown), 
                 method = "L-BFGS-B")
  }
  
  
  estimated_freqs <- fit$par
  hessian <- fit$hessian
  ses <- sqrt(diag(solve(hessian)))
  res <- data.frame(ancestry = colnames(prop_mat)[-c(1:n_known)], 
                    estimated_freq = estimated_freqs, 
                    low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                    high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
  
  return(res)
  
  
}
