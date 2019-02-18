#https://grunwaldlab.github.io/metacoder_documentation/example.html
#https://github.com/grunwaldlab/metacoder/issues/141  ## note our "otu_table" = "tax_table" in the example
#https://grunwaldlab.github.io/metacoder_documentation/workshop--05--plotting.html
#https://github.com/grunwaldlab/metacoder/issues/252

library(metacoder)
library(phyloseq)

setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")

psf_R<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psf_R.rds") # already rarified in Phyloseq
psf_Rp <- parse_phyloseq(psf_R)
print(psf_Rp) #pay attention to different tables and names

### Remove low count##
psf_Rp$data$otu_table <- zero_low_counts(psf_Rp, data = "otu_table", min_count = 2)

no_reads <- rowSums(psf_Rp$data$otu_table[, psf_Rp$data$sample_data$sample_id]) == 0
sum(no_reads)

psf_Rp <- filter_obs(psf_Rp, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(psf_Rp)

psf_Rp$data$tax_abund <- calc_taxon_abund(psf_Rp, "otu_table",
                                       cols = psf_Rp$data$sample_data$sample_id)

psf_Rp$data$tax_occ <- calc_n_samples(psf_Rp, "tax_abund", groups = psf_Rp$data$sample_data$Symptom_condition)

print(psf_Rp)

heat_tree(psf_Rp,
          node_size = rowMeans(cbind(Symptom, No_Symptom)),
          node_color =  Symptom - No_Symptom,
          node_label = taxon_names,
          tree_label = taxon_names,
          node_color_range = diverging_palette(),
          #node_color_range = c("black", "#FFFFFF"), #use for change main color
          #node_color_interval = c(-0.1, 0.1),
          #edge_color_interval = c(-0.1, 0.1),
          #edge_color_range = c("black", "#FFFFFF"),
          #node_label_max = 200,
          layout = "da", initial_layout = "re",
          title = "Symptom vs No-symptom",
          node_size_axis_label = "Average proportion",
          output_file = "metacoder_all_SvsNS.pdf")


## Get partial tree ##

psf_Rp %>%
  taxa::filter_taxa(taxon_names == "Ascomycota", subtaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  taxa::filter_taxa(taxon_ranks == "Family", supertaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = Symptom - No_Symptom,
            node_size_axis_label = "Average proportion",
            node_color_range = diverging_palette(),
            layout = "da", initial_layout = "re",
            output_file = "Asco_family_SvsNS.pdf")


####################
#### Proportion #### Doesn't really need here
####################
psf<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psf.rds") # already rarified in Phyloseq

psf_p <- parse_phyloseq(psf)
print(psf_p) #pay attention to different tables and names
psf_p$data$otu_table <- zero_low_counts(psf_p, data = "otu_table", min_count = 2)
no_reads <- rowSums(psf_p$data$otu_table[, psf_p$data$sample_data$sample_id]) == 0
sum(no_reads)
psf_p <- filter_obs(psf_p, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(psf_p)

psf_p$data$otu_table <- calc_obs_props(psf_p, "otu_table")


psf_p$data$tax_abund <- calc_taxon_abund(psf_p, "otu_table",
                                          cols = psf_p$data$sample_data$sample_id)

psf_p$data$tax_occ <- calc_n_samples(psf_p, "tax_abund", groups = psf_p$data$sample_data$Symptom_condition)

print(psf_p)

png("metacoder_proportion.png", width=4000, height=4000, res=400)
heat_tree(psf_p,
          node_size = rowMeans(cbind(Symptom, No_Symptom)),
          #node_size = n_obs,
          node_color =  Symptom - No_Symptom,
          node_label = taxon_names,
          tree_label = taxon_names,
          node_color_range = diverging_palette(),
          #node_color_range = c("black", "#FFFFFF"), #use for change main color
          #node_color_interval = c(-0.1, 0.1),
          #edge_color_interval = c(-0.1, 0.1),
          #edge_color_range = c("black", "#FFFFFF"),
          #node_label_max = 200,
          layout = "da", initial_layout = "re",
          title = "Symptom vs No-symptom",
          node_size_axis_label = "Average proportion"
          )
dev.off()
