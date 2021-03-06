###########################
##### Permanova table #####
###########################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(vegan)

adonis(otu_merged[,32:NCOL(otu_merged)] ~ STATUS/as.factor(Individual), method = "bray", data = otu_merged)
adonis(otu_merged[,32:NCOL(otu_merged)] ~ as.factor(HOUSE), method = "bray", data = otu_merged)
adonis(otu_merged[,32:NCOL(otu_merged)] ~ SUBSTATUS, method = "bray", data = otu_merged)

adonis(otu_merged[,32:NCOL(otu_merged)] ~ STATUS/as.factor(Individual) + as.factor(HOUSE) + SUBSTATUS, method = "bray", data = otu_merged)

adonis(otu_merged[,32:NCOL(otu_merged)] ~ as.factor(HOUSE) + STATUS/as.factor(Individual), method = "bray", data = otu_merged)

varpart(Y = otu_merged[,32:NCOL(otu_merged)], X = ~HOUSE, ~STATUS, ~as.factor(Individual), data = otu_merged)
