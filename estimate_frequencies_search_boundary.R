# Function to estimate ancestry-specific allele frequencies for a given variant.
# takes a vector of allele_counts allele counts for a given variant for n people 
# and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
# at the locus (could be genome-wide) for d ancestries. 
# this is a wrapper function that changes the frequency boundaries and 
# if frequencies for some ancestries are estimated at the set boundary conditions,
# compares the likelihood to the likelihood when the frequencies are set at the 
# boundaries of the parameter space (0 or 1, as needed).
estimate_frequencies_search_boundary <- function(allele_counts, prop_mat, confidence = 0.95, 
                                 frequency_boundary_grid = c(0.001, 0.01, 0.02, 0.05),
                                 use_smoothing_data = TRUE,
                                 chromosome_x = FALSE,
                                 sex = NULL, 
                                 male_label = "M", 
                                 mac_filter = 5){
  
  # prepare data.frame wor return null results if the function didn't ran
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
    return(return_val_not_run)
  }
  
  # start loop to find boundaries in which the analysis converges:
  converged <- FALSE
  ind_grid <- 1
  while (!converged & ind_grid <= length(frequency_boundary_grid)){
    low_freq_bound <- frequency_boundary_grid[ind_grid]
    high_freq_bound <- 1 - frequency_boundary_grid[ind_grid]
    ## optimize to estimate frequencies
    
    cur_res <- tryCatch({estimate_frequencies(allele_counts, prop_mat, confidence = confidence, 
                                    low_freq_bound = low_freq_bound, 
                                    high_freq_bound = high_freq_bound,
                                    use_smoothing_data = use_smoothing_data,
                                    chromosome_x = chromosome_x,
                                    sex = sex, 
                                    male_label = male_label, 
                                    mac_filter = mac_filter)},
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
  
  
  if (!converged) return(return_val_not_run)
  
  # check if convergence is right at the boundary condition
  low_bound_inds <- which(cur_res$res$estimated_freq == low_freq_bound)
  high_bound_inds <- which(cur_res$res$estimated_freq == high_freq_bound)
  
  if (length(low_bound_inds) > 0 | length(high_bound_inds) > 0){
  # if frequencies for some ancestries were estimated right at the lower bound:
  if (length(low_bound_inds) >0){
    known_freqs <- rep(0, length = length(low_bound_inds))
    names(known_freqs) <- res$ancestry[low_bound_inds]
  }
  
  # if frequencies for some ancestries were estimated right at the high bound:
  if (length(high_bound_inds) >0 ){
    known_freqs_high <- rep(1, length = length(high_bound_inds))
    names(known_freqs_high) <- res$ancestry[high_bound_inds]
    
    # should we augment a previous vector, and have a new known_freqs vector?
    if (length(low_bound_inds) > 0){
      known_freqs <- c(known_freqs, known_freqs_high)
    } else{
      known_freqs <- known_freqs_high
    }
  }
  
  # if some frequencies were estimated right on the provided boundaries, 
  # check if the likelihood is actually maximized at the boundary of the 
  # parameter space.
  if (length(known_freqs) > 0){
    boundary_res <- estimate_frequencies_w_known_freqs(allele_counts, prop_mat, 
                                                       known_freqs, confidence, 
                                                       low_freq_bound = frequency_boundary_grid[1], 
                                                       high_freq_bound = 1-frequency_boundary_grid[1],
                                                       use_smoothing_data = use_smoothing_data,
                                                       chromosome_x = chromosome_x,
                                                       sex = sex, 
                                                       male_label = male_label, 
                                                       mac_filter = mac_filter)
    
    if (boundary_res$nll < cur_res$nll) {
      # we need to return the new result, with the estimated 
      # frequency for the ancestry with boundary value set to that value
      return_val <- return_val_not_run
      return_val[match(boundary_res$res$ancestry, return_val$ancestry),] <- boundary_res$res
      return_val[match(names(known_freqs), return_val$ancestry),"estimated_freq"] <- known_freqs
      return(return_val)
    }
  }  
  
  return(cur_res$res)
}}





