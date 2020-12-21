
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
    
    # check that male counts are not higher than one (imputed data may have fractions)
    # if some dosages are higher than 1, divided all dosages/counts by 2. 
    max_male_dosage <- max(male_counts)
    if (max_male_dosage > 1){
      message("some male allele counts/dosages are larger than 1, dividing all counts by two...")
      male_counts <- male_counts/2
    }
    
    female_counts <- allele_counts[which(sex != male_label)]
    female_props <- prop_mat[which(sex != male_label),]
    
    # merge male and females:
    prop_mat <- rbind(male_props, female_props)
    allele_counts <- c(male_counts, female_counts)
    max_counts <- c(rep(1, length(male_counts)), rep(2, length(female_counts)))
  } else{  # not chromosome X
    
    # just add max_count
    max_counts <- rep(2, nrow(prop_mat))
    
  }
  
  return(list(allele_counts = allele_counts, 
              prop_mat = prop_mat, 
              max_counts = max_counts))
  
}
