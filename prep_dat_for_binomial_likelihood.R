
.prep_dat_for_binomial_likelihood <- function(allele_counts, prop_mat,  
                                            chromosome_x = FALSE,
                                            sex = NULL, 
                                            male_label = "M"){
  max_counts <- rep(NA, length(allele_counts))
  # If alleles are from chromosome X, need to be separately handled for males and females.
  if (chromosome_x){
    stopifnot(length(allele_counts) == length(sex))
    male_counts <- allele_counts[which(sex == male_label)]
    male_props <- prop_mat[which(sex == male_label),]
    
    # check that male counts are only zeros and ones, and not two.
    # if there are counts of two, turn them to one.
    inds_two <- which(male_counts == 2)
    if (length(inds_two) > 0){
      message("some male allele counts are equal to 2, setting them to 1...")
      male_counts[inds_two] <- 1
    }
    
    female_counts <- allele_counts[which(sex != male_label)]
    female_props <- prop_mat[which(sex != male_label),]
    
    # merge male and females:
    prop_mat <- rbind(male_props, female_props)
    allele_counts <- c(male_counts, females_counts)
    max_counts <- c(rep(1, nrow(male_counts)), rep(2, nrow(female_counts)))
  } else{  # not chromosome X
    
    # just add max_count
    max_counts <- rep(2, nrow(prop_mat))
    
  }
  
  return(list(allele_counts = allele_counts, 
              prop_mat = prop_mat, 
              max_counts = max_counts))
  
}
