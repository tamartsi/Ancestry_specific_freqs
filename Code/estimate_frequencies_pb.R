
# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 

### attempt to use the Poisson-Binomial distribution with LAFA ### 
### if it will work, we will merge it with the orignal function ###


estimate_frequencies_pb <- function(allele_counts, prop_mat, confidence = 0.95, 
                                 low_freq_bound = 0.001, high_freq_bound = 0.999,
                                 use_smoothing_data = FALSE,
                                 chromosome_x = FALSE,
                                 sex = NULL, 
                                 male_label = "M", 
                                 mac_filter = 5){
  stopifnot(length(allele_counts) == nrow(prop_mat))
  
  prop_mat <- as.matrix(prop_mat)
  
  # check for NAs, if there are observations with missing values, remove them.
  inds_na_alleles <- which(is.na(allele_counts))
  inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
  inds_na <- c(inds_na_alleles, inds_na_prop)
  if (length(inds_na) >0){
    message(paste(length(inds_na), "observations with missing values, removing them..."))
    allele_counts <- allele_counts[-inds_na]
    prop_mat <- prop_mat[-inds_na,]
    sex <- sex[-inds_na]
  }
  
  
  # check that we indeed need a separate prep function...
  prep_dat <- .prep_dat_for_poisbin_likelihood(allele_counts, prop_mat,
                                              chromosome_x = chromosome_x,
                                              sex = sex, 
                                              male_label = male_label)
  
  ancestry1 <- prep_dat$ancestry1
  ancestry2 <- prep_dat$ancestry2
  allele_counts <- prep_dat$allele_counts

  
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
  
  
  ## update to function that uses Poisson-Binomial 
  ## if this works we will merge it better and add "lafa" for a method
  ## as an argument.
  return(.estimate_frequencies_after_prep_pb(allele_counts, 
                                          prop_mat, 
                                          max_counts,
                                          low_freq_bound,
                                          high_freq_bound,
                                          confidence))
}


## separate function that is called after checks and preparations were done

# ancestry1 is a vector giving name of ancestry on one chromosomal copy 
# (no meaning for which copy) for participants. 
# ancestry1 is a vector giving name of ancestry on a second chromosomal copy
# if there is one (no meaning for which copy) for participants. 

.estimate_frequencies_after_prep_pb <- function(allele_counts,
                                             ancestry1,
                                             ancestry2 = NULL, # maybe at the end
                                          #   max_counts, # not needed
                                             low_freq_bound,
                                             high_freq_bound,
                                             confidence){
  
  ancestry_names <- unique(c(ancestry1, ancestry2))
  n_ancestry <- length(ancestry_names)
  ## compute the negative log likelihood function
  nll <- function(freqs){
    names(freqs) <- ancestry_names
    prob1 <- freqs[ancestry1]
    # prob1 and prob2 should be a function of freqs...
    nll_by_obs <- log(dpoisbinom_approx(allele_counts, prob1, prob2))
    return(-sum(nll_by_obs))		
  }
  
  
  ## optimize to estimate frequencies
  fit <- optim(par = rep(0.5, n_ancestry), fn = nll, hessian = TRUE, 
               lower = rep(low_freq_bound ,n_ancestry), 
               upper = rep(high_freq_bound, n_ancestry), 
               method = "L-BFGS-B")
  
  estimated_freqs <- fit$par
  hessian <- fit$hessian
  ses <- sqrt(diag(solve(hessian)))
  res <- data.frame(ancestry = ancestry_names, 
                    estimated_freq = estimated_freqs, 
                    low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                    high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
  
  # return the estimated frequencies, and the negative log likelihood. 
  return(list(res = res, nll = nll(res$estimated_freq)))
  
}


