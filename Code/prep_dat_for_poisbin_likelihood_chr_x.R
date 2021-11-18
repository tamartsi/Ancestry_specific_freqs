
# This function prepares data for using LAFA with 
# Poisson Binomial likelihood
# when the variant is on chromosome x

# instead of using "prop_mat", we now directly use vectors ancestry1, ancestry2.
# these vectors provide ancestry in one chromosomal copy, and in a second chromosomal
# copy. They are not aligned with allele counts (we don't know which ancestry 
# an allele is from).

# this function assumes that males only have ancestry1, and ancestry2 = NA for males.

.prep_dat_for_poisbin_likelihood_chr_x <- function(allele_counts, 
                                             ancestry1, 
                                             ancestry2 = NULL,
                                             sex, 
                                             male_label = "M"){
  
    stopifnot(length(allele_counts) == length(sex))
    male_inds <- which(sex == male_label)
  
    # check that male counts are not higher than one (imputed data may have fractions)
    # if some dosages are higher than 1, divided all dosages/counts by 2. 
    max_male_dosage <- max(allele_counts[male_inds])
    if (max_male_dosage > 1){
      message("some male allele counts/dosages are larger than 1, dividing all counts by two...")
      allele_counts[male_inds] <- allele_counts[male_inds]/2
    }  
    
    ## if all individuals are males, we don't need ancestry2 at all:
    if (length(male_inds) == length(allele_counts)){
      return(list(allele_counts = allele_counts, 
                  ancestry1 = ancestry1, 
                  ancestry2 = NULL))
    }
    
    ## if there are some females and some males, set ancestry2 of males to NA
    ancestry2[male_inds] <- NA  
    
    return(list(allele_counts = allele_counts, 
                ancestry1 = ancestry1, 
                ancestry2 = ancestry2))

}
