
.prep_dat_for_binary_likelihood <- function(allele_counts, prop_mat,  
                                            chromosome_x = FALSE,
                                            sex = NULL, 
                                            male_label = "M"){
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
    
    decomposed_fem_alleles <- sapply( female_counts,decompose_two_alleles_one_person)
    decomposed_fem_alleles <- as.numeric(matrix(decomposed_fem_alleles)	)
    female_props_double <- female_props[rep(1:nrow(female_props), each = 2),]
    
    # merge male and females:
    prop_mat_double <- rbind(male_props, female_props_double)
    decomposed_alleles <- c(male_counts, decomposed_fem_alleles)
  } else{  # not chromosome X
    
    # turn each allele dosage into a two rows, each reporting one allele
    # with the same proportion ancestries.
    decomposed_alleles <- sapply( allele_counts,decompose_two_alleles_one_person)
    decomposed_alleles <- as.numeric(matrix(decomposed_alleles)	)
    prop_mat_double <- prop_mat[rep(1:nrow(prop_mat), each = 2),]
    
  }
  
  return(list(decomposed_alleles = decomposed_alleles, 
              prop_mat_double = prop_mat_double))
  
}
