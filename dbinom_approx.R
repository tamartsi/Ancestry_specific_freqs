
# function that returns binomial probability mass function, 
# with linear interpolation to the case where we have dosage data.

dbinom_approx <- function(x, size, prob){
  arg2 <- prob^(x)*(1-prob)^(size-x)
  
  # now compute arg1
  floor_nums <- floor(x) 
  ceiling_nums <- ceiling(x)
  
  floor_vals <- ifelse(floor_nums == 1, 2, 1)
  ceiling_vals <- ifelse(ceiling_nums == 1, 2, 1)
  
  arg1 <- floor_vals + (ceiling_vals - floor_vals)*(x-floor_nums)
  
  return(arg1*arg2)
}