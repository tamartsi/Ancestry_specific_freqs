require(data.table)

choose_obs <- function(x){
  return(x[sample(1:length(x), size =1)])
}

# simulate data by first generating global ancestry proportion, then computing 
# the expected value of allele probability per person, and sampling from this
# probability
simulate_data <- function(eff_n_anc, maf_anc ,lai_length =0.05, names_anc = NULL, 
                           error_rate = 0){
  stopifnot(floor(sum(eff_n_anc)) - sum(eff_n_anc) == 0 , length(maf_anc) == length(eff_n_anc))
  stopifnot( all(floor(eff_n_anc/lai_length) - eff_n_anc/lai_length  ==0))
  # LAI are of fixed length of 5% of the genome.
  
  if (is.null(names_anc)){
    names_anc <- paste0("anc_", 1:length(eff_n_anc))
  }
  n <- sum(eff_n_anc)
  n_lai <- n/lai_length
  lai_df <- data.table(lai_num = 1:n_lai, 
                       lai_anc = do.call(c, mapply(rep, names_anc, eff_n_anc/lai_length)), 
                       person_id = paste0("person_", sample(rep(1:n, each = 1/lai_length ), size = n/lai_length)))
  setkeyv(lai_df, "person_id")
  
  # compute proportion global ancestries for each person:  
  prop_df <- lai_df[, length(lai_num)*lai_length, by = c("person_id", "lai_anc")]
  setnames(prop_df, "V1", "prop_anc")
  # reformat 
  prop_df_wide <- dcast(prop_df, person_id ~ lai_anc, value.var = "prop_anc")
  for (anc in names_anc){
    na_inds <- which(is.na(prop_df_wide[[anc]]))
    prop_df_wide[[anc]][na_inds] <- 0
  } 
  
  
  # focus on two LAIs (one from each chromosomal copy).
  # first, choose at random the two LAIs.
  lai_copy_1_df <- lai_df[, choose_obs(lai_num), by = c("person_id")]
  setnames(lai_copy_1_df, "V1", "lai_num")
  lai_copy_1_df <- merge(lai_copy_1_df, lai_df[, c("lai_num", "lai_anc")], by = "lai_num")
  
  lai_df_rm_first <- lai_df[-which(lai_df$lai_num %in% lai_copy_1_df$lai_num)]
  
  lai_copy_2_df <- lai_df_rm_first[, choose_obs(lai_num), by = c("person_id")]
  setnames(lai_copy_2_df, "V1", "lai_num")
  lai_copy_2_df <- merge(lai_copy_2_df, lai_df[, c("lai_num", "lai_anc")], by = "lai_num")
  
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
  
  
  # In case we want to induce errors, generate another dataset with errors:
  if (error_rate >0){
    sample_inds_err_copy1 <- sample(1:nrow(dat), size= floor(nrow(dat)*error_rate))
    sample_inds_err_copy2 <- sample(1:nrow(dat), size= floor(nrow(dat)*error_rate))
    
    dat$lai_copy1_anc_err <- dat$lai_copy1_anc
    dat$lai_copy2_anc_err <- dat$lai_copy2_anc
    for (i in 1:length(sample_inds_err_copy1)){
      cur_ind <- sample_inds_err_copy1[i]
      cur_anc <- dat$lai_copy1_anc_err[cur_ind]
      dat$lai_copy1_anc_err[cur_ind] <- sample(setdiff( names_anc, cur_anc), 1)
      
      cur_anc <- dat$lai_copy2_anc_err[cur_ind]
      dat$lai_copy2_anc_err[cur_ind] <- sample(setdiff( names_anc, cur_anc), 1)
    }
    
    for (anc in names_anc){
      dat[[paste0(anc, "_count_err")]] <- as.numeric(dat$lai_copy1_anc_err == anc) + 
        as.numeric(dat$lai_copy2_anc_err == anc)
    }
    
  }
  

  
  ## now return three datasets:
  # prop_anc: proportion ancestries
  # allele_count: genotype allele count
  # lai_count: lai allele counts, when the genotype from allele_count is in the lai. 
  
  # if a dataset with error was requested, return this too

  prop_anc <- as.data.frame(dat[,c("person_id", names_anc), with = FALSE])
  allele_count <- as.data.frame(dat[,c("person_id", "allele_count"), with = FALSE])
  lai_count <- as.data.frame(dat[,c("person_id", paste0(names_anc, "_count")), with = FALSE])
  
  if (error_rate == 0){
    return(list(prop_anc = prop_anc, 
                allele_count = allele_count, 
                lai_count = lai_count))
  } else{
    lai_count_err <- as.data.frame(dat[,c("person_id", paste0(names_anc, "_count_err")), with = FALSE])
    return(list(prop_anc = prop_anc, 
                allele_count = allele_count, 
                lai_count = lai_count, 
                lai_count_err = lai_count_err))
  }
  
  
}

# check
eff_n_anc <- c(9, 11, 80)
maf_anc <- c(0.5, 0.3, 0.2)

dat <- simulate_data(eff_n_anc , maf_anc)

dat <- simulate_data(eff_n_anc , maf_anc, error_rate = 0.05)
