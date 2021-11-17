require(poisbinom)

# function that returns poisson binomial probability mass function, 
# with linear interpolation to the case where we have dosage data.

# the function receives 
# x: a count of alleles
# prob1: allele frequency of one ancestral allele
# prob2: allele frequency of second ancestral allele. 
# prob2 is NULL if there is only one allele (e.g. X-chromosome in a male-only sample)



dpoisbinom_approx <- function(x, prob1, prob2 = NULL){
  
  # check length match
  stopifnot(length(x) == length(prob1))
  
  # checks to account for numerical error
  if (sum(prob1[!is.na(prob1)] > 1 + 1e-10 | prob1[!is.na(prob1)] < -1e-10) > 0)  
    stop("prob1 higher than 1 or lower than 0")
  
  prob1 <- ifelse(prob1 > 1, 1, prob1)
  
  dat <- cbind(x, prob1)
  
  # same checks for prob2 -- only if not NULL (usually it won't be)
  if (!is.null(prob2)){
    stopifnot(length(x) == length(prob2))
    
    if (sum(prob2[!is.na(prob2)] > 1 + 1e-10 | prob2[!is.na(prob2)] < -1e-10) > 0)  
      stop("prob1 higher than 1 or lower than 0")
    
    prob2 <- ifelse(prob2 < 0 , 0, prob2)
    
    dat <- cbind(dat, prob2)
  }
  
  if (is.null(prob2)) stopifnot(all(x>= 0 & x<=1))
  
  # now apply dpoisbinom_approx_1 on the data
  est_probs <- apply(dat, 1, function(x){
                    dpoisbinom_approx_1(x[1], x[2], x[3])
                      })
  return(est_probs)
  
}


# for one person, x, prob1, prob2 are scalars (prob2 may be NULL)
# assume values are checked and clean.
dpoisbinom_approx_1 <- function(x, prob1, prob2 = NULL){
  
  prob_vec <- na.omit(c(prob1, prob2))
  
  # now compute arg1
  floor_x <- floor(x) 
  ceiling_x <- ceiling(x)
  
  floor_val <- dpoisbinom(floor_x, prob_vec)
  ceiling_val <- dpoisbinom(ceiling_x, prob_vec)
  
  floor_val + (ceiling_val - floor_val)*(x-floor_x)
  
}
