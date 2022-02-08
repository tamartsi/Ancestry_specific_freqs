require(data.table)


# this function updates simulat_data_keep_phase by setting a specific proportion
# of ancestry-heterozygote individuals.


# simulate data by first generating global ancestry proportion, then computing 
# the expected value of allele probability per person, and sampling from this
# probability.
# unlike in the function simulate_data, we here keep the phasing information 
# (i.e. to know which ancestry is which allele). This will be used to study 
# the increase in estimation efficiency when knowing the ancestry of each allele.
simulate_data_keep_phase_enrich_het_ancestry <- function(eff_n_anc, 
                                                         maf_anc,
                                                         prop_het, 
                                                         lai_length =0.05, 
                                                         names_anc = NULL, 
                                                         error_rate = 0){
  stopifnot(floor(sum(eff_n_anc)) - sum(eff_n_anc) == 0 , length(maf_anc) == length(eff_n_anc))
  stopifnot( all(floor(eff_n_anc/lai_length) - eff_n_anc/lai_length  ==0))
  
  # LAI are of fixed length of 5% (lai_length*100) of the genome.
  num_lai_in_genome <- 1/lai_length
  
  if (is.null(names_anc)){
    names_anc <- paste0("anc_", 1:length(eff_n_anc))
  }
  n <- sum(eff_n_anc)
  names(eff_n_anc) <- names_anc
  
  # simulate the first LAI with a fixed proportion of ancestry-heterozygote individuals.
  heteros <- sample(1:n, size = round(n*prop_het))
  non_heteros <- setdiff(1:n, heteros)
  # LAI 1: choose ancestries at random
  lai_1 <- sample(names_anc, size = n, replace = TRUE)
  
  # LAI 2 (which will be the equivalent to LAI 1 on the second chromosomal copy)
  # initialize, then will out values making sure we have hetero and homozygotes as needed.
  lai_2 <- rep(NA, length = n)
  lai_2[non_heteros] <- lai_1[non_heteros]
  for (i in heteros){
    lai_2[i] <- sample(setdiff(names_anc, lai_1[i]), 1)
  }
  
  lai_12_df <- data.table(lai_num = c(rep("lai_1", n), rep("lai_2", n)),
                          lai_anc = c(lai_1, lai_2), 
                          person_id = rep(paste0("person_", 1:n) , 2))
  
  # compute how many LAIs we need from each ancestry
  tot_lai <- eff_n_anc/lai_length
  need_lai <- rep(NA, length(eff_n_anc))
  for (i in 1: length(need_lai)){
    need_lai[i] <- tot_lai[i] - sum(lai_1 == names(eff_n_anc)[i]) - sum(lai_2 == names(eff_n_anc)[i])
  }
  
  # total remaining  LAIs to simulate to have the correct effN for each ancestry:
  need_lai_vec <- do.call(c, mapply(rep, names_anc, need_lai))
 
  lai_df <- data.table(lai_num = paste0("lai_", rep(3:num_lai_in_genome, n)), 
                       lai_anc = need_lai_vec, 
                       person_id = rep(paste0("person_", 1:n), num_lai_in_genome -2))
  
  lai_df <- rbind(lai_df, lai_12_df)
  
  setkeyv(lai_df, "person_id")
  
  
  # compute proportion global ancestries for each person:  
  prop_df <- lai_df[, length(lai_num)/num_lai_in_genome, by = c("person_id", "lai_anc")]
  setnames(prop_df, "V1", "prop_anc")
  # reformat 
  prop_df_wide <- dcast(prop_df, person_id ~ lai_anc, value.var = "prop_anc")
  for (anc in names_anc){
    na_inds <- which(is.na(prop_df_wide[[anc]]))
    prop_df_wide[[anc]][na_inds] <- 0
  } 
  
  
  # focus on LAI 1 and LAI 2 (corresponding to each other each chromosomal copy).
  lai_copy_1_df <- lai_df[lai_num == "lai_1"]
  lai_copy_2_df <- lai_df[lai_num == "lai_2"]
  
  
  setnames(lai_copy_1_df, "lai_anc", "lai_copy1_anc")
  setnames(lai_copy_2_df, "lai_anc", "lai_copy2_anc")
  lai_df_per_person <- merge(lai_copy_1_df[, c("person_id", "lai_copy1_anc"), with = FALSE],
                             lai_copy_2_df[, c("person_id", "lai_copy2_anc"), with = FALSE])
  


  ## merge data to a single data.table:
  dat <- merge(lai_df_per_person, prop_df_wide, by = "person_id")
  
  
  ## now add frequencies for each chromosomal copy, then generate variants, and sum:
  table_anc_freq <- data.table(lai_copy1_anc = names_anc, 
                               lai_copy2_anc = names_anc, 
                               lai_copy1_freq = maf_anc, 
                               lai_copy2_freq = maf_anc)
  
  dat <- merge(dat, table_anc_freq[,c("lai_copy1_anc", "lai_copy1_freq"), with = FALSE], 
               by = "lai_copy1_anc", keep = "all")
  
  dat <- merge(dat, table_anc_freq[,c("lai_copy2_anc", "lai_copy2_freq"), with = FALSE], 
               by = "lai_copy2_anc", keep = "all")
  
  dat[, allele_copy1 := rbinom(1, size = 1, prob = lai_copy1_freq) , by = "person_id"]
  dat[, allele_copy2 := rbinom(1, size = 1, prob = lai_copy2_freq), by = "person_id" ]
  
  dat[,allele_count := allele_copy1 + allele_copy2]
  
  ## add lai counts by ancestry, to study a use case where we now the local 
  ## ancestry count in the interval overlapping the genotype:
  for (anc in names_anc){
    dat[[paste0(anc, "_count")]] <- as.numeric(dat$lai_copy1_anc == anc) + 
      as.numeric(dat$lai_copy2_anc == anc)
  }
  

  global_prop_anc <- as.data.frame(dat[,c("person_id", names_anc), with = FALSE])
  allele_by_copy <- as.data.frame(dat[,c("person_id", "allele_copy1", "allele_copy2"), with = FALSE])
  ancestry_by_copy <- as.data.frame(dat[,c("person_id", "lai_copy1_anc", "lai_copy2_anc")])
  
  allele_count <- as.data.frame(dat[,c("person_id", "allele_count"), with = FALSE])
  lai_count <- as.data.frame(dat[,c("person_id", paste0(names_anc, "_count")), with = FALSE])
  
  # compute local proportion ancestry:
  local_prop_anc <- lai_count
  colnames(local_prop_anc) <- sub("_count", "", colnames(lai_count))
  local_prop_anc[,names_anc] <- local_prop_anc[,names_anc]/rowSums(local_prop_anc[,names_anc])

  return(list(global_prop_anc = global_prop_anc, 
              local_prop_anc = local_prop_anc,
              allele_by_copy = allele_by_copy,
              ancestry_by_copy = ancestry_by_copy,
              allele_count = allele_count, 
              lai_count = lai_count))

}



