
##
# Function to estimate ancestry-specific allele frequencies for a given variant,
# when some (one or more) ancestral allele frequencies are known and provided.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# known_freq need to be a named vector, specifying known frequencies. The names of elements 
# in known_freq are some of the column names of prop_mat. 



LAFA_w_known_freqs <- function(allele_counts, 
                               ancestry1, 
                               ancestry2 = NULL,
                               known_freqs, 
                               confidence = 0.95, 
                               low_freq_bound = 0.001, high_freq_bound = 0.999,
                               use_smoothing_data = FALSE,
                               chromosome_x = FALSE,
                               sex = NULL, 
                               male_label = "M",
                               mac_filter = 5,
                               ancestry_names = NULL){
  
  stopifnot(length(allele_counts) == length(ancestry1))
  if (is.null(ancestry_names)) ancestry_names <- unique(c(ancestry1, ancestry2))
  
  stopifnot(all(is.element(names(known_freqs), ancestry_names)))
  stopifnot(length(known_freqs) > 0 & length(known_freqs) < length(ancestry_names))
  

  
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
    smoothing_data <- .generate_smoothing_observations_lafa(setdiff(ancestry_names, names(known_freqs)))
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    ancestry1 <- c(ancestry1, smoothing_data$simulated_ancestry1)
    
    if (!is.null(ancestry2)){
      ancestry2 <- c(ancestry2, rep(NA, length(smoothing_data$simulated_ancestry1)))
    }
  }
  

  
  out <- .LAFA_w_known_freq_after_prep(allele_counts, 
                                       ancestry1, 
                                       ancestry2 = ancestry2,
                                       known_freqs, 
                                       n_unknown,
                                       n_known,
                                       low_freq_bound, 
                                       high_freq_bound,
                                       confidence)
  out$res <- out$res[match(ancestry_names, out$res$ancestry),]
  out$res$ancestry <- ancestry_names
  out$res$estimated_freq[match(names(known_freqs), out$res$ancestry)] <- known_freqs
  return(out)
  
}



.LAFA_w_known_freq_after_prep <- function(allele_counts, 
                                         ancestry1, 
                                         ancestry2 = NULL,
                                         known_freqs, 
                                         n_unknown,
                                         n_known,
                                         low_freq_bound, 
                                         high_freq_bound,
                                         confidence){
  
  ancestry_names <- na.omit(unique(c(ancestry1, ancestry2)))
  unknown_ancestry_names <- setdiff(ancestry_names, names(known_freqs))
  n_ancestry <- length(ancestry_names)
  n_unknown <- length(unknown_ancestry_names)
 
  ## compute the negative log likelihood function
  
  nll <- function(unknown_freqs){
    names(unknown_freqs) <- unknown_ancestry_names
    freqs <- c(known_freqs, unknown_freqs)
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
  res <- data.frame(ancestry = unknown_ancestry_names, 
                    estimated_freq = estimated_freqs, 
                    low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                    high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
  
  # return the estimated frequencies, and the negative log likelihood. 
  return(list(res = res, nll = nll(res$estimated_freq)))
  
}


