
# Function to estimate ancestry-specific allele frequencies for a given variant.
# now we assumes that we have local ancestry, and the ancestry calls are phased 
# relative to the allele calls.
# for each of n people, we have the ancestry in their two copies, and alleles in their two copies.
# for now, we allele_counts has two columns, one for each chromosomal copy; 
# similarly alleles_ancestry. This can be easily redesigned.
# we assume that the order of people in the same in allele_counts and in alleles_ancestry.
# no checks this time...
estimate_frequencies_phased <- function(allele_counts, allele_ancestry, 
                                        confidence = 0.95, 
                                 low_freq_bound = 0.001, high_freq_bound = 0.999,
                                 use_smoothing_data = FALSE,
                                 mac_filter = 5, 
                                 allele_copy_label = "copy"){
  stopifnot(nrow(allele_counts) == nrow(allele_ancestry))
 
  copy1_col_allele_counts <- grep(paste0(allele_copy_label, "1"), colnames(allele_counts), value = TRUE)
  copy2_col_allele_counts <- grep(paste0(allele_copy_label, "2"), colnames(allele_counts), value = TRUE)
  
  copy1_col_allele_ancestry <- grep(paste0(allele_copy_label, "1"), colnames(allele_ancestry), value = TRUE)
  copy2_col_allele_ancestry <- grep(paste0(allele_copy_label, "2"), colnames(allele_ancestry), value = TRUE)
  
  allele_counts_formatted <- c(allele_counts[,copy1_col_allele_counts], allele_counts[,copy2_col_allele_counts])
  allele_ancestry_formatted <- c(allele_ancestry[,copy1_col_allele_ancestry], allele_ancestry[,copy2_col_allele_ancestry])
  
  ancestry_names <- unique(allele_ancestry_formatted)
  prop_mat <- matrix(0, ncol = length(ancestry_names), nrow = length(allele_ancestry_formatted))
  colnames(prop_mat) <- ancestry_names
  for (anc in ancestry_names){
    prop_mat[which(allele_ancestry_formatted == anc),anc] <- 1
  }
    
 
  max_counts <- rep(1, length(allele_counts_formatted))
  
  # check if the number of minor alleles is higher than the mac_filter,
  # stop if MAC is too low.
  stopifnot(min(sum(allele_counts_formatted), sum(max_counts) - sum(allele_counts_formatted)) > mac_filter)
  
  # add made-up data to avoid boundaries of the frequency parameter space   
  if (use_smoothing_data){
    smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
    prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat)
    allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
    max_counts <- c(max_counts, rep(1, nrow(smoothing_data$simulated_prop_mat)))
  }
  
  
  
  return(.estimate_frequencies_after_prep(allele_counts_formatted, 
                                          prop_mat, 
                                          max_counts,
                                          low_freq_bound,
                                          high_freq_bound,
                                          confidence))
}

