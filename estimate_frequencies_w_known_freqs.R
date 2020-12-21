## UPDATE: using Binomial distribution instead of Bernoulli.
##
# Function to estimate ancestry-specific allele frequencies for a given variant,
# when some (one or more) ancestral allele frequencies are known and provided.
# takes a vector of allele_countss allele counts for a given variant for n people 
# and an n x d matrix matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# known_freq need to be a named vector, specifying known frequencies. The names of elements 
# in known_freq are some of the column names of prop_mat. 
estimate_frequencies_w_known_freqs <- function(allele_counts, prop_mat, known_freqs, confidence = 0.95, 
                                            low_freq_bound = 0.001, high_freq_bound = 0.999,
                                            use_smoothing_data = FALSE,
                                            chromosome_x = FALSE,
                                            sex = NULL, 
                                            male_label = "M",
                                            mac_filter = 5){
  
  stopifnot(length(allele_counts) == nrow(prop_mat))
  stopifnot(all(is.element(names(known_freqs), colnames(prop_mat))))
  stopifnot(length(known_freqs) > 0 & length(known_freqs) < ncol(prop_mat))
  
  prop_mat <- as.matrix(prop_mat)
  
  # check for NAs, if there are observations with missing values, remove them.
  inds_na_alleles <- which(is.na(allele_counts))
  inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
  inds_na <- c(inds_na_alleles, inds_na_prop)
  if (length(inds_na) >0){
    message(paste(length(inds_na), "observations with missing values, removing them..."))
    allele_counts <- allele_counts[-inds_na]
    prop_mat <- prop_mat[-inds_na,]
  }
  
  # re-order the column so that the ones with known frequencies are at the beginning:
  n_known <- length(known_freqs)
  n_unknown <- ncol(prop_mat) - n_known
  inds_known <- which(is.element(colnames(prop_mat), names(known_freqs)))
  inds_unknown <- setdiff(1:ncol(prop_mat), inds_known)
  prop_mat <- prop_mat[,c(inds_known, inds_unknown)]
  
  prep_dat <- .prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
                                              chromosome_x = chromosome_x,
                                              sex = sex, 
                                              male_label = male_label)
  
  prop_mat <- prep_dat$prop_mat
  allele_counts <- prep_dat$allele_counts
  max_counts <- prep_dat$max_counts
  
  # check if the number of minor alleles is higher than the mac_filter,
  # stop if MAC is too low.
  stopifnot(min(sum(allele_counts), sum(max_counts) - sum(allele_counts)) > mac_filter)

  # add made-up data to avoid boundaries of the frequency parameter space 
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
    # add only the rows corresponding to the ancestries with unknown frequencies
    inds_keep <- (n_known*2 + 1):((n_known + n_unknown)*2)
    prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat[inds_keep,])
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count[inds_keep])
    max_counts <- c(max_counts, rep(1, length(inds_keep)))
  }
  
  
  
  known_sums <- prop_mat[,names(known_freqs), drop = FALSE] %*% known_freqs
  
  
  
  ## compute the negative log likelihood function
  nll <- function(unknown_freqs){
    allele_probs <- as.numeric(known_sums + prop_mat[,c(n_known+1):ncol(prop_mat), drop = FALSE] %*% unknown_freqs)
    nll_by_obs <- log(dbinom_approx(allele_counts, max_counts, allele_probs))
    return(-sum(nll_by_obs))		
  }
  
  ## optimize to obtain allele frequency estimates
  
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






