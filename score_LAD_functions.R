#####
## R function to compute LAD score 
####

## This workflow calculates ancestry probabilities for individuals based on genetic differences to a set of proxy haplotypes. 
## The method uses a weighting scheme where individuals genetically closer to a proxy haplotype receive higher probability contributions for the proxy’s region.

library(dplyr)
library(tidyr)

## INPUT:
## table_diff_SAS_meta : 
#This is a table where: Each row represents an individual in the study. 
#First columns represent meta data (SampleID, LP genotype, Region of origin).
#Then columns represent proxy haplotypes (e.g., I11456_1, I8726_2). They are store in 'proxy_haplotype' of them
#Values are integer differences between the individual’s haplotype and the proxy haplotype.

## metadata_proxy:
#Table with information of the name of proxy haplotype and their region or origin

###
## Convert data to right format
###

##Converts wide format (one column per proxy haplotype) into long format.
long_table_diff_SAS_meta <- table_diff_SAS_meta %>%
  pivot_longer(
    cols = 4:(length(proxy_haplotype)+3),          # all the value columns (e.g., I11456_1, I8726_2, etc.)
    names_to = "Proxy_Ind_hap",              # Name of the proxy haplotype.
    values_to = "Nb_diff"               # Number of differences between the individual and the proxy haplotype.
  )

##Adding information about the 'ancestral source' of the proxy (Region_proxy) from 'metadata_proxy' table
long_table_diff_SAS_meta=merge(long_table_diff_SAS_meta,metadata_proxy,by.x="Proxy_Ind_hap",by.y="SampleID",suffixes = c("","_proxy"))

###
## Score computation
###
lambda=10 #weight
ancestry_probs <- long_table_diff_SAS_meta %>%
  mutate(weight = 1 / (1 + lambda * Nb_diff)) %>%
  group_by(SampleID, LCT_genotype, Region_proxy) %>%
  # average across proxies for each region (avoids bias from regions with more proxies)
  summarise(region_mean = mean(weight), .groups = "drop") %>%
  group_by(SampleID, LCT_genotype) %>%
  # normalize so probabilities sum to 1 per individual
  mutate(prob = region_mean / sum(region_mean)) %>%
  ungroup()

## Filtering with ancestry thresold
ancestry_probs80 <- filter(ancestry_probs,prob>0.80)
