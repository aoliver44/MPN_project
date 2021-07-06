##### GFLASSO ######

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(tidyverse)
library(corrplot)
library(pheatmap)
library(gflasso)
library(janitor)
library(viridis)

cytokines <- read.csv("Cytokines.txt", sep = "\t")
cytokines_mpn <- cytokines %>% filter(., Status == "MPN", Sample != "124", Sample != "65") 
mpn_w_cytokine <- as.vector(unique(cytokines_mpn$Sample))


mpn_otu <- subset(otu_merged, otu_merged$Individual %in% mpn_w_cytokine)
rownames(mpn_otu) <- NULL
mpn_otu <- mpn_otu %>% 
  select(., Individual, 32:NCOL(.)) %>% 
  group_by(., Individual) %>%
  summarise_all(funs(mean)) 
mpn_otu <- mpn_otu %>% `row.names<-`(., NULL) %>% column_to_rownames(var = "Individual") %>% clean_names()
mpn_otu <- mpn_otu[,(colSums(mpn_otu) > 2)]

cytokines_mpn <- cytokines_mpn[cytokines_mpn$Sample %in% rownames(mpn_otu), ]
cytokines_mpn <- cytokines_mpn %>% select(., 1, 3:NCOL(cytokines_mpn)) %>% 
  clean_names() 
rownames(cytokines_mpn) <- cytokines_mpn$sample
cytokines_mpn$sample <- NULL

DS <- cor(cytokines_mpn, method = "spearman")
corrplot(DS)

#CV <- readRDS(file = "CV_gflasso.rds")
#system.time(CV <- cv_gflasso(X = scale(mpn_otu), Y = scale(cytokines_mpn), R = DS, nCores = 4, 
#                             additionalOpts = list(delta_conv = 1e-5, iter_max = 1e5)))

# cv_plot_gflasso(CV)
# 
#gfMod <- gflasso(X = scale(mpn_otu), Y = scale(cytokines_mpn), R = DS, opts = list(lambda = CV$optimal$lambda,
#                                                                    gamma = CV$optimal$gamma, 
#                                                                    delta_conv = 1e-5,
#                                                                    iter_max = 1e5))

gfMod <- readRDS(file = "gfMod.rds")
colnames(gfMod$B) <- colnames(cytokines_mpn)
Lasso_data <- gfMod$B[abs(rowSums(gfMod$B)) > 0.3, ]
rownames(Lasso_data) <- sapply(strsplit(rownames(Lasso_data), '_g_'), getElement, 2)
pheatmap(Lasso_data, show_rownames = T)
heatmap(Lasso_data, scale = "row", col=viridis(40), margins = c(4,15))


# check raw data
cyto_otu <- merge(cytokines_mpn, mpn_otu, by.x = "row.names", by.y = "row.names")

ggplot(data = cyto_otu, aes(y = dialister, x = phascolarctobacterium)) +
  geom_point()

ggplot(data = cyto_otu, aes(y = tn_fa, x = parabacteroides)) +
         geom_point()

ggplot(data = cyto_otu, aes(y = tn_fa, x = cyto_otu$k_bacteria_p_actinobacteria_c_coriobacteriia_o_coriobacteriales_f_coriobacteriaceae_g_collinsella_s_stercoris)) +
  geom_point()      

# Corr with cytokines
library(corrr)
cyto_otu %>% select(., tn_fa, 16:86) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(tn_fa) %>%
  filter(abs(tn_fa) > 0.1) %>% 
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  ggplot(aes(x = rowname, y = tn_fa)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (R) with Fiber\n (Spearman)") +
  xlab("Species") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

library(rstatix)
tnf_cor <- cyto_otu %>% select(., tn_fa, 16:86) %>%
  cor_test(method = "spearman", vars = "tn_fa")
tnf_cor <- tnf_cor %>% filter(., var1 == "tn_fa") 
tnf_cor1 <- as.data.frame(p.adjust(tnf_cor$p, method = "fdr"))
tnf_cor <- cbind(tnf_cor, tnf_cor1)
