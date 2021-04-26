# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# this is a wrapper function that changes the frequency boundaries and 
# if frequencies for some ancestries are estimated at the set boundary conditions,
# compares the likelihood to the likelihood when the frequencies are set at the 
# boundaries of the parameter space (0 or 1, as needed).
estimate_frequencies_dynamic_boundary <- function(allele_counts, prop_mat, confidence = 0.95, 
                                 frequency_boundary_grid = c(0.001, 0.01, 0.02, 0.05),
                                 use_smoothing_data = TRUE,
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
  
  # prepare data.frame for returning null results if the function doesn't run
  # or does not converge
  return_val_not_run <- data.frame(ancestry = colnames(prop_mat),
                                   estimated_freq = NA,
                                   low_CI = NA,
                                   high_CI = NA)
  
   # check that MAC is higher than mac_filter
  prep_dat <- .prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
                                                chromosome_x = chromosome_x,
                                                sex = sex, 
                                                male_label = male_label)
  
  prop_mat <- prep_dat$prop_mat
  allele_counts <- prep_dat$allele_counts
  max_counts <- prep_dat$max_counts
  

  if (min(sum(allele_counts), sum(max_counts) - sum(allele_counts)) <= mac_filter){
    return(list(res = return_val_not_run, 
                nll = NA,
                boundary = NA, 
                flag = "MAC filter"))
  }
  
  # add made-up data to avoid boundaries of the frequency parameter space   
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
    prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat)
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    max_counts <- c(max_counts, rep(1, nrow(smoothing_data$simulated_prop_mat)))
  }
  
  
  # start loop to find boundaries in which the analysis converges:
  converged <- FALSE
  ind_grid <- 1
  while (!converged & ind_grid <= length(frequency_boundary_grid)){
    low_freq_bound <- frequency_boundary_grid[ind_grid]
    high_freq_bound <- 1 - frequency_boundary_grid[ind_grid]
    ## optimize to estimate frequencies
    
    cur_res <- tryCatch({.estimate_frequencies_after_prep(allele_counts = allele_counts, 
                                                          prop_mat = prop_mat, 
                                                          max_counts = max_counts,
                                                          freqs = freqs,
                                                          low_freq_bound = low_freq_bound,
                                                          high_freq_bound = high_freq_bound,
                                                          confidence = confidence)},
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
  # prepare a vector called set_freqs with in which the ancestries
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
  if (length(set_freqs) == ncol(prop_mat)){
    allele_probs <- prop_mat[,names(set_freqs), drop = FALSE] %*% set_freqs
    nll_by_obs <- log(dbinom_approx(allele_counts, max_counts, allele_probs))
    nll <- -sum(nll_by_obs)
    
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
    tryCatch({.estimate_frequencies_w_known_freq_after_prep(allele_counts = allele_counts, 
                                                  prop_mat = prop_mat, 
                                                  max_counts = max_counts,
                                                  known_freqs= set_freqs, 
                                                  n_unknown = ncol(prop_mat) - length(set_freqs),
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





