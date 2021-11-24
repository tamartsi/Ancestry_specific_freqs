# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# this is a wrapper function that changes the frequency boundaries and 
# if frequencies for some ancestries are estimated at the set boundary conditions,
# compares the likelihood to the likelihood when the frequencies are set at the 
# boundaries of the parameter space (0 or 1, as needed).

### LAFA update #### (need to update documentation above; check/test)

LAFA_dynamic_boundary <- function(allele_counts, 
                                  ancestry1, 
                                  ancestry2 = NULL,
                                  confidence = 0.95, 
                                 frequency_boundary_grid = c(0.001, 0.01, 0.02, 0.05),
                                 use_smoothing_data = TRUE,
                                 chromosome_x = FALSE,
                                 sex = NULL, 
                                 male_label = "M", 
                                 mac_filter = 5, 
                                 ancestry_names = NULL){
 
  
  stopifnot(length(allele_counts) == length(ancestry1))
  if (is.null(ancestry_names)) ancestry_names <- unique(c(ancestry1, ancestry2))
  
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
  
  
  
  
  
  
  # add made-up data to avoid boundaries of the frequency parameter space   
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations_lafa(unique(c(ancestry1, ancestry2)))
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    ancestry1 <- c(ancestry1, smoothing_data$simulated_ancestry1)
    
    if (!is.null(ancestry2)){
      ancestry2 <- c(ancestry2, rep(NA, length(smoothing_data$simulated_ancestry1)))
    }
  }
  
  # prepare data.frame for returning null results if the function doesn't run
  # or does not converge
  return_val_not_run <- data.frame(ancestry = ancestry_names,
                                   estimated_freq = NA,
                                   low_CI = NA,
                                   high_CI = NA)
  
 
 
  
  # start loop to find boundaries in which the analysis converges:
  converged <- FALSE
  ind_grid <- 1
  while (!converged & ind_grid <= length(frequency_boundary_grid)){
    low_freq_bound <- frequency_boundary_grid[ind_grid]
    high_freq_bound <- 1 - frequency_boundary_grid[ind_grid]
    ## optimize to estimate frequencies
    
    cur_res <- tryCatch({.estimate_frequencies_after_prep_pb(allele_counts, 
                                                             ancestry1,
                                                             ancestry2,
                                                             low_freq_bound,
                                                             high_freq_bound,
                                                             confidence)},
                        error = function(error){
                          "not converged"
                        }
                      )
    if (!identical(cur_res, "not converged")){
      converged <- TRUE
    } else{
      ind_grid <- ind_grid + 1
    }
  }
  
  
  if (!converged) {
    return(list(res = return_val_not_run, 
                nll = NA,
                boundary = NA, 
                flag = "Not converged"))
  } else{ # if converged, prepare potential return value
    cur_res$res <- cur_res$res[match(ancestry_names, cur_res$res$ancestry),]
    cur_res$res$ancestry <- ancestry_names
    potential_return_val_converged <- list(res = cur_res$res,
                                           nll = cur_res$nll,
                                           boundary = low_freq_bound,
                                           flag = "converged")
  }
  
  # check if convergence is right at the boundary condition for some ancestries
  low_bound_inds <- which(cur_res$res$estimated_freq == low_freq_bound)
  high_bound_inds <- which(cur_res$res$estimated_freq == high_freq_bound)
  
  # if nothing was estimated at the boundary, return the computed frequencies.
  if (length(low_bound_inds) + length(high_bound_inds) == 0){
    return(potential_return_val_converged)
  }
  
  # otherwise, we continue to identify which frequencies were 
  # estimated at the boundaries and...
  # prepare a vector called set_freqs in which the ancestries
  # where frequencies were estimated in the boundaries are set to the 
  # boundaries of the parameter space. Start with low and then high boundary.
    
  # if frequencies for some ancestries were estimated right at the lower bound:
  if (length(low_bound_inds) >0){
      set_freqs <- rep(0, length = length(low_bound_inds))
      names(set_freqs) <- cur_res$res$ancestry[low_bound_inds]
    }
  
  # if frequencies for some ancestries were estimated right at the high bound:
  if (length(high_bound_inds) >0 ){
      set_freqs_high <- rep(1, length = length(high_bound_inds))
      names(set_freqs_high) <- cur_res$res$ancestry[high_bound_inds]
      
      # should we augment a previous vector, and have a new set_freqs vector?
      if (length(low_bound_inds) > 0){
        set_freqs <- c(set_freqs, set_freqs_high)
      } else{
        set_freqs <- set_freqs_high
      }
    } # finished preparing a vector of set frequencies
 
  # first: if all frequencies were estimated at the boundary:
  if (length(set_freqs) == length(ancestry_names)){
    prob1 <- set_freqs[ancestry1]
    if (!is.null(ancestry2)){
      prob2 <- set_freqs[ancestry2]
    }  else{
      prob2 <- NULL
    }
    like_by_obs <- dpoisbinom_approx(allele_counts, prob1, prob2)

    
    # check if we have an impossible value (probability zero)
    if (sum(like_by_obs == 0) > 0) return(potential_return_val_converged)
    
    nll <- -sum(log(like_by_obs))
    
    if (nll < potential_return_val_converged$nll){
      return_val <- return_val_not_run
      return_val[match(names(set_freqs), return_val$ancestry),"estimated_freq"] <- set_freqs
      return(list(res = return_val,
                  nll = nll,
                  boundary = potential_return_val_converged$boundary,
                  flag = "converged and evaluated monomorphic"))
    } else{
      return(potential_return_val_converged)
    }
  }
  
  # at this point: some frequencies but not all were estimated at the boundary
  boundary_res <- 
    tryCatch({.LAFA_w_known_freq_after_prep(allele_counts, 
                                            ancestry1, 
                                            ancestry2 = NULL,
                                            known_freqs= set_freqs, 
                                            n_unknown = length(ancestry_names) - length(set_freqs),
                                            n_known = length(set_freqs),
                                            low_freq_bound = low_freq_bound, 
                                            high_freq_bound = high_freq_bound,
                                            confidence = confidence)},
             error = function(error){
               "not converged"
             }
    )
    # if this did not converge...
    if (identical(boundary_res, "not converged")){
     return(potential_return_val_converged)
    }
  
    # if it did converge, continue...
  
    if (boundary_res$nll < cur_res$nll) {
      # we need to return the new result, with the estimated 
      # frequency for the ancestry with boundary value set to that value
      return_val <- return_val_not_run
      return_val[match(boundary_res$res$ancestry, return_val$ancestry),] <- boundary_res$res
      return_val[match(names(set_freqs), return_val$ancestry),"estimated_freq"] <- set_freqs
      return(list(res = return_val,
                  nll = boundary_res$nll,
                  boundary = low_freq_bound,
                  flag = "converged some freqs evaluated monomorphic")
                  )
    } else{ 
      return(potential_return_val_converged)
    } 

}





