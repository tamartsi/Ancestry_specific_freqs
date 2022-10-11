# About this repository 
The repository provides code for "AFA: Ancestry-specific allele frequency estimation in admixed populations: The Hispanic Community Health Study/Study of Latinos", a method described in [this paper](https://www.sciencedirect.com/science/article/pii/S2666247722000124). "AFA" computes ancestry-specific allele frequencies in admixed populations using maximum likelihood estimation. The method can use either local proportion ancestries, "LAFA", or global ancestry proportions "GAFA". The repository also provides a dataset of U.S. Hispanic/Latino ancestry-specific allele frequencies and their confidence intervals estimated based on the Hispanic Community Health Study/Study of Latinos (HCHS/SOL) using "GAFA" and "LAFA". HCHS/SOL is an admixed population with three predominant ancestries: Amerindian, European, and African. The estimates are based on HCHS/SOL imputed data. 

# Contact information: 
Tamar Sofer:  tsofer@bwh.harvard.edu
Einat Granot-Hershkovitz:  egranot-hershkovitz@bwh.harvard.edu

# Code
This folder provides functions to compute ancestry-specific allele frequencies, as well as code for simulations, and some tests/examples. 

# HCHS_SOL_GAFA_estimates
This folder provides compressed files for 23 chromosomes with U.S. Hispanic/Latino ancestry-specific allele frequencies and their confidence intervals estimated based on the HCHS/SOL using "GAFA". "GAFA" method uses global ancestry estimates: estimated proportion of genomes inherited from each of the ancestral population used. 

# HCHS_SOL_LAFA_estimates
This folder provides compressed files for 23 chromosomes with U.S. Hispanic/Latino ancestry-specific allele frequencies and their confidence intervals estimated based on the HCHS/SOL using "LAFA". "LAFA" method uses local ancestry interval estimates which were previously inferred using the RFMix software. Overall, 15,500 local ancestry intervals are dispersed throughout the genome and are used as proportions of genetic ancestry for variants located in the intervals. 
Local ancestry intervals (LAIs)
Three-way LAI (Amerindian, African, and European) were previously inferred in 12,793 HCHS/SOL individuals using the RFMix software with a reference panel derived from the combination of the Human Genome Diversity Project (HGDP) and the 1000 Genome Project (using the GRCh37 assembly) representing the relevant ancestral populations20. Overall, 15,500 are LAIs dispersed throughout the genome (14,815 LAI in autosomal chromosomes),
# Column annotation of the compressed 23 chromosomal files providing ancestry-specific allele frequencies:

Each file has a few columns, as follows: 
- CHR:	Chromosome number
- POS:	Chromosomal position (hg38)
- allele_a:	Reference allele
- allele_b:	Alternative allele
- Africa_estimated_freq:	Estimated African frequency 
- Africa_low_CI:	Low point of the 95% CI of the African frequency
- Africa_high_CI:	High point of the 95% CI of the African frequency
- Europe_estimated_freq:	Estimated European frequency
- Europe_low_CI:	Low point of the 95% CI of the European frequency 
- Europe_high_CI:	High point of the 95% CI of the European frequency
- America_estimated_freq:	Estimated Amerindian frequency 
- America_low_CI:	Low point of the 95% CI of the Amerindian frequency
- America_high_CI:	High point of the 95% CI of the Amerindian frequency
- boundary:	Boundary condition used when the function converges. 
- imputed: Imputed genotype
- R2:	Imputation quality

# Acknowledgements
We thank the staff and participants of the Hispanic Community Health Study/Study of Latinos for their important contributions.
