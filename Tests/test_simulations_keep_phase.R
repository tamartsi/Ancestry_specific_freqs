# check
eff_n_anc <- c(9, 11, 80)
maf_anc <- c(0.5, 0.3, 0.2)

dat <- simulate_data_keep_phase(eff_n_anc , maf_anc)


# use phasing of ancestry relative to variant alleles when estimating frequencies:
estimate_frequencies_phased(allele_counts = dat$allele_by_copy, 
                            allele_ancestry = dat$ancestry_by_copy)

# without using the phasing of the ancestry relative to the variant alleles,
# using global ancestries:
estimate_frequencies(allele_counts = dat$allele_count$allele_count, 
                      prop_mat = dat$global_prop_anc[, c("anc_1", "anc_2", "anc_3")])


# without using the phasing of the ancestry relative to the variant alleles,
# using local ancestries:
estimate_frequencies(allele_counts = dat$allele_count$allele_count, 
                     prop_mat = dat$local_prop_anc[, c("anc_1", "anc_2", "anc_3")])
