
# first test the function providing the density mass function for poisson bionomial 
# distribution for given frequencies under LAFA

# for one person:
# note this function assumes that values are checked (so it doesn't check that 
# x cannot be larger than the number of allele copies, 
# i.e. of probabilities provided)
dpoisbinom_approx_1(x = 1, prob1 = 0.5)
dpoisbinom_approx_1(x = 1, prob1 = 0.5, prob2 = 1)
dpoisbinom_approx_1(x = 0, prob1 = 0.5, prob2 = 1)
dpoisbinom_approx_1(x = 2, prob1 = 0.5, prob2 = 1)
dpoisbinom_approx_1(x = 1.8, prob1 = 0.5, prob2 = 1)

dpoisbinom_approx_1(x = 1, prob1 = 0.5, prob2 = NA)

# for dpoisbinom_approx that uses vectorized data (n people at once)
# this one does some check

## should have an error:
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        3), 
                      prob1 = c(0.2, 0.2, 0.6, 0.6,
                                0.6, 0.6, 0.2, 0.2, 
                                0.2), 
                      prob2 = c(0.2, 0.2, 0.2, 0.2, 
                                0.6, 0.2, 0.2, 0.6,
                                0.6))
# should be fine
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, 0.2, 0.2, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            0.6))

# suppose prob2 is NULL, but some x values are 2:
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2))


# suppose prob2 is NULL, x values are fine
dpoisbinom_approx(x = c(1, 0,0,0,
                        1,0,0,0,
                        1), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2))


# some prob2 values are NAs:
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, NA, NA, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            0.6))
# some prob2 values are NAs but when x = 2:
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, 0.2, 0.2, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            NA))

# implausible probability values:
dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, 0.2, 0.2, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            1.1))


# numerical problem -- should be fixed by the function
res1 <- dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, 0.2, 0.2, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            -1e-10))

res2 <- dpoisbinom_approx(x = c(1, 0,0,0,
                        2,0,0,0,
                        2), 
                  prob1 = c(0.2, 0.2, 0.6, 0.6,
                            0.6, 0.6, 0.2, 0.2, 
                            0.2), 
                  prob2 = c(0.2, 0.2, 0.2, 0.2, 
                            0.6, 0.2, 0.2, 0.6,
                            0))
identical(res1, res2)


.estimate_frequencies_after_prep_pb(allele_counts = c(1, 0,0,0,
                                                      2,0,0,0,
                                                      2),
                                    ancestry1 = c("a1", "a1", "a2", "a2", 
                                                  "a1", "a2", "a2", "a2", 
                                                  "a2"), 
                                    ancestry2 = c("a1", "a1", "a2", "a2",
                                                  "a1", "a1", "a2", "a2", 
                                                  "a1"),
                                    low_freq_bound  = 0.01,
                                    high_freq_bound = 0.99,
                                    confidence = 0.95)
# it worked! 


# this needs to change counts of 2 in males to 1, 
# and set ancestry2 for males to NA.
.prep_dat_for_poisbin_likelihood_chr_x(allele_counts = c(1, 0,0,0,
                                                         2,0,0,0,
                                                         2),
                                       ancestry1 = c("a1", "a1", "a2", "a2", 
                                                     "a1", "a2", "a2", "a2", 
                                                     "a2"),
                                       ancestry2 = c("a1", "a1", "a2", "a2",
                                                     "a1", "a1", "a2", "a2", 
                                                     "a1"),
                                       sex = c("F", "F", "F", "F",
                                               "M", "F", "F", "M",
                                               "M"), 
                                       male_label = "M")


# all are males, ancestry2 is set to NULL, all counts are divided by half.
.prep_dat_for_poisbin_likelihood_chr_x(allele_counts = c(1, 0,0,0,
                                                         2,0,0,0,
                                                         2),
                                       ancestry1 = c("a1", "a1", "a2", "a2", 
                                                     "a1", "a2", "a2", "a2", 
                                                     "a2"),
                                       ancestry2 = c("a1", "a1", "a2", "a2",
                                                     "a1", "a1", "a2", "a2", 
                                                     "a1"),
                                       sex = c("M", "M", "M", "M",
                                               "M", "M", "M", "M",
                                               "M"), 
                                       male_label = "M")


estimate_frequencies_pb(allele_counts = c(1, 2,0,0,
                                          2,0,0,0,
                                          2),
                        ancestry1 = c("a1", "a1", "a2", "a2", 
                                      "a1", "a2", "a2", "a2", 
                                      "a2"),
                        ancestry2 = c("a1", "a1", "a2", "a2",
                                      "a1", "a1", "a2", "a2", 
                                      "a1"))


estimate_frequencies_pb(allele_counts = c(1, 2,0,0,
                                          2,0,0,0,
                                          2, 1, 1,2),
                        ancestry1 = c("a1", "a1", "a2", "a2", 
                                      "a1", "a2", "a2", "a2", 
                                      "a2", "a1", "a1", "a2"),
                        ancestry2 = c("a1", "a1", "a2", "a2",
                                      "a1", "a1", "a2", "a2", 
                                      "a1", "a1", "a2", "a2"), 
                        use_smoothing_data = TRUE)
