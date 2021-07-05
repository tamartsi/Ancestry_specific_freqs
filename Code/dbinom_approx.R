
# function that returns binomial probability mass function, 
# with linear interpolation to the case where we have dosage data.

dbinom_approx <- function(x, size, prob){
  
  # add checks to account for numerical error
  if (sum(prob > 1 + 1e-10 | prob < -1e-10) > 0)  
    stop("probability higher than 1 or lower than 0")
  prob <- ifelse(prob > 1, 1, prob)
  prob <- ifelse(prob <0 , 0, prob)
  
  arg2 <- prob^(x)*(1-prob)^(size-x)
  
  # now compute arg1
  floor_nums <- floor(x) 
  ceiling_nums <- ceiling(x)
  
  floor_vals <- ifelse(floor_nums == 1, 2, 1)
  ceiling_vals <- ifelse(ceiling_nums == 1, 2, 1)
  
  arg1 <- floor_vals + (ceiling_vals - floor_vals)*(x-floor_nums)
  
  return(arg1*arg2)
}
