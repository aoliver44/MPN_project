##########################
### GENERATE BASIC ENV ###
##########################

#path_to_data <- ""
#setwd(path_to_data)

library(EcolUtils)

# Import OTU file and rarefy
qiime_raw <- t(read.csv("feature-table-MPN.tsv", row.names = 1, check.names = FALSE, sep = "\t"))
qiime_clean <- qiime_raw
barplot(sort(rowSums(qiime_clean)), ylim = c(0, max(rowSums(qiime_clean))), 
        xlim = c(0, NROW(qiime_clean)), col = "Blue", main = "Reads per Sample") 
sort(rowSums(qiime_clean))

# rarefy to 2200 sequences
set.seed(seed = 999)
rare_perm_otu <- rrarefy.perm(qiime_clean, sample = 10000, n = 100, round.out = T)
alpha_rare <- rare_perm_otu[rowSums(rare_perm_otu) >= 10000-(10000*.1), colSums(rare_perm_otu) >= 1]
barplot(sort(rowSums(alpha_rare)), ylim = c(0, max(rowSums(rare_perm_otu))), 
        xlim = c(0,NROW(rare_perm_otu)), col = "Blue", main = "Reads after rarefaction")
sort(rowSums(alpha_rare))

# add taxonomic information
library(tidyverse)
library(stringr)
library(zoo)
otu_with_taxa <- as.data.frame(alpha_rare)
taxonomy <- read.csv("taxonomy.tsv", sep = "\t")
taxonomy$Taxon <- gsub("[", "", taxonomy$Taxon, fixed = T)
taxonomy$Taxon <- gsub("]", "", taxonomy$Taxon, fixed = T)

# what this next part does is attempt to bring the last known phylgenetic group to unknown sections
taxonomy <- taxonomy %>% separate(., col = Taxon, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop")
taxonomy <- data.frame(lapply(taxonomy, function(x) {gsub("[a-z]__", "", x)}))
taxonomy <- taxonomy %>% mutate_all(na_if,"")
taxonomy <- data.frame(t(apply(taxonomy,1,function(x) na.locf(x))))

taxonomy$taxonomy <- paste0("k__", taxonomy$L1,";","p__", taxonomy$L2,";","c__", taxonomy$L3,";","o__", taxonomy$L4,";","f__", taxonomy$L5,";","g__", taxonomy$L6,";","s__", taxonomy$L7)
names(otu_with_taxa) <- taxonomy$taxonomy[match(names(otu_with_taxa), taxonomy$Feature.ID)]
names(otu_with_taxa) <- make.unique(colnames(otu_with_taxa))
#write.csv(otu_with_taxa, file = "rarefied_otu_table.txt", sep = "\t")

# import metadata
metadata <- read.csv("metadata_2.14.19.txt", check.names = F, sep = "\t")

# merge OTU with taxa and metadata
otu_merged <- merge(metadata, otu_with_taxa, by.x = "Sample", by.y = "row.names") 
otu_merged$SAMPLE <- rownames(otu_merged)
otu_merged$SAMPLE <- NULL
#unload all packages and clean up before next analysis:
rm(rare_perm_otu)
rm(taxonomy)
rm(qiime_raw)
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

