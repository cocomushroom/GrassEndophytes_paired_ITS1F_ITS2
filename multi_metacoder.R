library(metacoder)
library(phyloseq)
setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")

psf_R_en<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psf_R_en.rds") # already rarified in Phyloseq

psf_Renp <- parse_phyloseq(psf_R_en)
print(psf_Renp) #pay attention to different tables and names

psf_Renp$data$otu_table <- zero_low_counts(psf_Renp, data = "otu_table", min_count = 2)

no_reads <- rowSums(psf_Renp$data$otu_table[, psf_Renp$data$sample_data$sample_id]) == 0
sum(no_reads)

psf_Renp <- filter_obs(psf_Renp, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(psf_Renp)

psf_Renp$data$tax_abund <- calc_taxon_abund(psf_Renp, "otu_table",
                                          cols = psf_Renp$data$sample_data$sample_id)

psf_Renp$data$tax_occ <- calc_n_samples(psf_Renp, "tax_abund", groups = psf_Renp$data$sample_data$Plant_species)
print(psf_Renp)


# Calculate difference between treatments
psf_Renp$data$diff_table <- compare_groups(psf_Renp, data = "otu_table",
                                    cols = psf_Renp$data$sample_data$sample_id,
                                    groups = psf_Renp$data$sample_data$Plant_species)

# Plot results (might take a few minutes)
png("metacoder_multi.png", width=6000, height=6000, res=400)
heat_tree_matrix(psf_Renp,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 #edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions")

dev.off()

