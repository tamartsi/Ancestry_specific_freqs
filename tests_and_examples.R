
prop_mat <- matrix(c(0.1, 0.7, 0.2, 0.2, 0.2, 0.6, 0.4, 0.2, 0.4, 0.2, 0.3, 0.5, 0.6, 0.2, 0.2, 1, 0, 0, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 8, byrow = TRUE)
colnames(prop_mat) <- c("eur", "afr", "amr")
rowSums(prop_mat) # all are 1. 
true_freqs <- c(eur =0.2, afr = 0.8, amr = 0.2)
probs <- prop_mat %*% true_freqs
allele_count <- rbinom(length(probs), 2, probs )
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.05, high_freq_bound = 0.95)

# this sometimes fails: 
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.0001, high_freq_bound = 0.9999) 
# with smoothing, this works:
estimate_frequencies(allele_count, prop_mat,low_freq_bound = 0.0001, high_freq_bound = 0.9999, use_smoothing_data = TRUE)


estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.1))

# with smoothing: 
estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.2), 
                                                       use_smoothing_data = TRUE)


estimate_frequencies_w_known_freqs(allele_count, prop_mat, known_freqs = c(eur=0.1, afr = 0.8))
