#https://grunwaldlab.github.io/metacoder_documentation/example.html
#https://github.com/grunwaldlab/metacoder/issues/141
library(metacoder)
library(phyloseq)


setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")

psf_R<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psf_R.rds") # already rarified in Phyloseq

psf_Rp <- parse_phyloseq(psf_R)
print(psf_Rp) #pay attention to different tables and names

psf_Rp$data$otu_table <- zero_low_counts(psf_Rp, data = "otu_table", min_count = 2)

no_reads <- rowSums(psf_Rp$data$otu_table[, psf_Rp$data$sample_data$sample_id]) == 0
sum(no_reads)

psf_Rp <- filter_obs(psf_Rp, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(psf_Rp)

psf_Rp$data$tax_abund <- calc_taxon_abund(psf_Rp, "otu_table",
                                       cols = psf_Rp$data$sample_data$sample_id)

psf_Rp$data$tax_occ <- calc_n_samples(psf_Rp, "tax_abund", groups = psf_Rp$data$sample_data$Symptom_condition)

print(psf_Rp)
png("metacoder.png", width=3000, height=3500,res=250)
heat_tree(psf_Rp,
          node_size = n_obs,
          node_color =  Symptom - No_Symptom,
          node_label = taxon_names,
          tree_label = taxon_names)
dev.off()




