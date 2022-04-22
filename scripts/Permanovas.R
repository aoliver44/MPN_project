###########################
##### Permanova table #####
###########################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(vegan)
library(tidyverse)
adonis(otu_merged[,32:NCOL(otu_merged)] ~ STATUS/as.factor(Individual), method = "bray", data = otu_merged) 

## HOUSE
otu_merged_tmp <- otu_merged %>% filter(., HOUSE != "")
otu_merged_tmp$HOUSE <- as.factor(otu_merged_tmp$HOUSE)
otu_merged_tmp$HOUSE <- droplevels(otu_merged_tmp$HOUSE)
levels(otu_merged_tmp$HOUSE)

adonis(otu_merged_tmp[,32:NCOL(otu_merged_tmp)] ~ STATUS/as.factor(Individual) + as.factor(HOUSE) , method = "bray", data = otu_merged_tmp)

## SUBSTATUS
otu_merged_tmp <- otu_merged %>% filter(., SUBSTATUS != "")
otu_merged_tmp$SUBSTATUS <- droplevels(otu_merged_tmp$SUBSTATUS)
levels(otu_merged_tmp$SUBSTATUS)

adonis(otu_merged_tmp[,32:NCOL(otu_merged_tmp)] ~ SUBSTATUS/as.factor(Individual), method = "bray", data = otu_merged_tmp)
