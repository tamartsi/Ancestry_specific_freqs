
# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 

### attempt to use the Poisson-Binomial distribution with LAFA ### 
### if it will work, we will merge it with the original function ###
### or create separate functions, one for GAFA, one for LAFA ###

# despite not being traditional notation to set default null values, 
# doing it here due to potential chromosome X variants. Most variants are not, though.

estimate_frequencies_pb <- function(allele_counts, 
                                    ancestry1, 
                                    ancestry2 = NULL, 
                                    confidence = 0.95, 
                                 low_freq_bound = 0.001, high_freq_bound = 0.999,
                                 use_smoothing_data = FALSE,
                                 chromosome_x = FALSE,
                                 sex = NULL, 
                                 male_label = "M", 
                                 mac_filter = 5){
  stopifnot(length(allele_counts) == length(ancestry1))
  
  
  # check for NAs, if there are observations with missing values, remove them.
  inds_na_alleles <- which(is.na(allele_counts))
  inds_na_ancestry1 <- which(is.na(ancestry1))
  
  inds_na <- c(inds_na_alleles, inds_na_ancestry1)
  
  if (!chromosome_x){
    inds_na_ancestry2 <- which(is.na(ancestry2))
    inds_na <- c(inds_na, inds_na_ancestry2)
  }
  
  
  if (length(inds_na) >0){
    message(paste(length(inds_na), "observations with missing values, removing them..."))
    allele_counts <- allele_counts[-inds_na]
    ancestry1 <- ancestry1[-inds_na]
    
    if (!chromosome_x){
      ancestry2 <- ancestry2[-inds_na]
    }
    
  }
  
  if (chromosome_x){
    stopifnot(length(sex) == length(allele_counts))
    prep_dat <- .prep_dat_for_poisbin_likelihood_chr_x(allele_counts, 
                                                 ancestry1,
                                                 ancestry2,
                                                 sex = sex, 
                                                 male_label = male_label)
    
    ancestry1 <- prep_dat$ancestry1
    ancestry2 <- prep_dat$ancestry2
    allele_counts <- prep_dat$allele_counts
    
  }
  
  

  
  # check if the number of minor alleles is higher than the mac_filter,
  # stop if MAC is too low.
  if (!chromosome_x){
    max_alleles <- 2*length(allele_counts)  
  } else{ # chromosome x
    n_male <- sum(sex == male_label)
    n_all <- length(allele_counts)
    max_alleles <- n_male + 2*(n_all - n_male)
  }
  stopifnot(min(sum(allele_counts), 
                max_alleles - sum(allele_counts)) > mac_filter)    
  
  
  # add made-up data to avoid boundaries of the frequency parameter space   
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations_lafa(unique(c(ancestry1, ancestry2)))
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    ancestry1 <- c(ancestry1, smoothing_data$simulated_ancestry1)
    
    if (!is.null(ancestry2)){
      ancestry2 <- c(ancestry2, rep(NA, length(smoothing_data$simulated_ancestry1)))
    }
  }
  
  
  ## update to function that uses Poisson-Binomial 
  ## if this works we will merge it better and add "lafa" for a method
  ## as an argument.
  return(.estimate_frequencies_after_prep_pb(allele_counts, 
                                          ancestry1,
                                          ancestry2,
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
                                             ancestry2 = NULL,
                                             low_freq_bound,
                                             high_freq_bound,
                                             confidence){
  
  ancestry_names <- na.omit(unique(c(ancestry1, ancestry2)))
  n_ancestry <- length(ancestry_names)
  ## compute the negative log likelihood function
  nll <- function(freqs){
    names(freqs) <- ancestry_names
    prob1 <- freqs[ancestry1]
    if (!is.null(ancestry2)){
      prob2 <- freqs[ancestry2]
    }  else{
      prob2 <- NULL
    }
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


