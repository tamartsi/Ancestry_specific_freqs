
# Take an allele count (values 0, 1, or 2) and return a vector of allele counts
# in two chromosomes (values 0 or 1)

decompose_two_alleles_one_person <- function(allele_count){
  
  # check that only valid allele counts are present
  stopifnot(is.element(allele_count, c(0,1,2)))
  
  if (allele_count == 0){
    allele_count_two <- c(0,0)
  }
  if (allele_count == 1){
    allele_count_two <- c(1,0)
  }
  if (allele_count == 2){
    allele_count_two <- c(1,1)
  }
  
  return(c(decomposed_allele_count =  allele_count_two))
}
