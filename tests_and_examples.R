
prop_mat <- matrix(c(0.1, 0.7, 0.2, 0.2, 0.2, 0.6, 0.4, 0.2, 0.4, 0.2, 0.3, 0.5, 0.6, 0.2, 0.2, 1, 0, 0, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 8, byrow = TRUE)
colnames(prop_mat) <- c("eur", "afr", "amr")
rowSums(prop_mat) # all are 1. 
true_freqs <- c(eur =0.2, afr = 0.8, amr = 0.2)
probs <- prop_mat %*% true_freqs
allele_count <- rbinom(length(probs), 2, probs )

# check internal function: 

.prep_dat_for_binomial_likelihood(allele_count, prop_mat)


estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.05, high_freq_bound = 0.95)
estimate_frequencies_search_boundary(allele_count, prop_mat)

true_freqs <- c(eur =0.2, afr = 0.8, amr = 0)
probs <- prop_mat %*% true_freqs
allele_count <- rbinom(length(probs), 2, probs )
estimate_frequencies_search_boundary(allele_count, prop_mat)

# this sometimes fails: 
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.0001, high_freq_bound = 0.9999) 
# with smoothing, this works:
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.0001, high_freq_bound = 0.9999, use_smoothing_data = TRUE)

# one frequency is 0:
true_freqs <- c(eur =0.2, afr = 0.8, amr = 0)
probs <- prop_mat %*% true_freqs
allele_count <- rbinom(length(probs), 2, probs )
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.0001, high_freq_bound = 0.9999, use_smoothing_data = TRUE)
estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(amr=0), mac_filter = 1)


estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.1))

# with missing values: 
prop_mat_na <- prop_mat
prop_mat_na[3,2] <- NA
estimate_frequencies_w_known_freqs(allele_count, prop_mat_na, known_freqs = c(eur=0.1))
estimate_frequencies(allele_count, prop_mat_na)

# missing values also of allele_count:
allele_count_na <- allele_count
allele_count_na[5] <- NA
estimate_frequencies_w_known_freqs(allele_count_na, prop_mat_na, known_freqs = c(eur=0.1))
estimate_frequencies(allele_count_na, prop_mat_na)


# with smoothing: 
estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.2), 
                                                       use_smoothing_data = TRUE)


estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.1, afr = 0.8))


# check X chromosome: 
prop_mat <- matrix(c(0.1, 0.7, 0.2, 0.2, 0.2, 0.6, 0.4, 0.2, 0.4, 0.2, 0.3, 0.5, 0.6, 0.2, 0.2, 1, 0, 0, 0, 1, 0, 0.1, 0.8, 0.1, 
                     0.4, 0.4, 0.2, 0.5, 0.3, 0.2), nrow = 10, byrow = TRUE)
colnames(prop_mat) <- c("eur", "afr", "amr")
rowSums(prop_mat) # all are 1. 
true_freqs <- c(eur =0.7, afr = 0.8, amr = 0.2)
probs <- prop_mat %*% true_freqs
allele_count <- rbinom(length(probs), 2, probs )
sex <- c("f", "m", "f", "f", "f", "m", "f", "f", "f", "m")
estimate_frequencies(allele_count, prop_mat,
                     low_freq_bound = 0.05, 
                     high_freq_bound = 0.95, 
                     chromosome_x = TRUE,
                    sex = sex, 
                    male_label = "m")
# if we don't provide sex...
estimate_frequencies(allele_count, prop_mat,
                     low_freq_bound = 0.05, 
                     high_freq_bound = 0.95, 
                     chromosome_x = TRUE,
                     male_label = "m")
# error! (as should be)