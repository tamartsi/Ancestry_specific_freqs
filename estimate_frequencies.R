## UPDATE: using Binomial distribution instead of Bernoulli.
##
# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
estimate_frequencies <- function(allele_counts, prop_mat, confidence = 0.95, 
                                 low_freq_bound = 0.001, high_freq_bound = 0.999,
                                 use_smoothing_data = FALSE,
                                 chromosome_x = FALSE,
                                 sex = NULL, 
                                 male_label = "M", 
                                 mac_filter = 5){
  stopifnot(length(allele_counts) == nrow(prop_mat))
  
  prop_mat <- as.matrix(prop_mat)
  
  # check for NAs, if there are observations with missging values, remove them.
  inds_na_alleles <- which(is.na(allele_counts))
  inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
  inds_na <- c(inds_na_alleles, inds_na_prop)
  if (length(inds_na) >0){
    message(paste(length(inds_na), "observations with missing values, removing them..."))
    allele_counts <- allele_counts[-inds_na]
    prop_mat <- prop_mat[-inds_na,]
  }
  
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
    prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat)
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    max_counts <- c(max_counts, rep(1, nrow(smoothing_data$simulated_prop_mat)))
  }
  
  
  ## compute the negative log likelihood function
  nll <- function(freqs){
    allele_probs <- as.numeric(prop_mat %*% freqs)
    nll_by_obs <- log(dbinom_approx(allele_counts, max_counts, allele_probs))
    return(-sum(nll_by_obs))		
  }
  
  
  ## optimize to estimate frequencies
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





