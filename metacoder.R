library(metacoder)
library(phyloseq)


setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")
psf<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psf.rds")

psf_p <- parse_phyloseq(psf)

png("metacoder.png", width=3000, height=3500,res=250)
heat_tree(psf_p,
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names)
dev.off()

