###############################################
########### Starting phyloseq##################
setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")
library("ape")
library("vegan")
library("tidyr")
library("metacoder")
library(phyloseq)
library("ggplot2")

## Load taxa produced by dada2##
tax<-readRDS(file = "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/tax_final.rds")
## Load seqtab.nochim produced by dada2##
seqtab.nochim<-readRDS(file = "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/seqtabnochim.rds")

## Not a QIIME2 format
meta <-read.delim("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/metadata_ITS_update.txt", row.names = 1)


## Load table for mycotoxin## ZEAR_4_SULF, Ergine and Avenacein_Y_Q1 are only presense/absense
toxin_value<-data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/value.txt"))
toxin_meta<- data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/metadata.txt"))

toxin_value<- apply(toxin_value, 2, function(x){gsub("BLOQ",0.5,x)}) ## Replace cells with peak on HPLC but below detection limit to 0.5
toxin_value[is.na(toxin_value)] <- 0
toxin_value[toxin_value==""]<-0
toxin_value<- data.frame(MycoID=toxin_value[,1], apply(toxin_value[,-1], 2, as.numeric))

toxin_quality<-toxin_value
toxin_quality[toxin_quality>0]<-1
toxin_quality$MycoID<-toxin_value$MycoID

## Merge with Mycotoxin#


##### Start Phyloseq#####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(tax))
ps
# Keep fungi only
psf <- subset_taxa(ps, Kingdom=="k__Fungi")

# check read distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(psf), TRUE), sorted = 1:ntaxa(psf), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(psf), 
                                                        TRUE), sorted = 1:nsamples(psf), type = "Samples"))
png("read_distribution.png", width=2000, height=2500,res=250)
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()

nreads = sort(sample_sums(psf)) # if need the lowest read number for a sample to determine the rarification number
nreads


#### Create Rarified table ### 
set.seed(1111)
psf_R = rarefy_even_depth(psf, sample.size = nreads[[1]])

## to check if the rarefaction did work
#par(mfrow = c(1, 2))
title = "Sum of reads for each sample, physeq_R"
plot(sort(sample_sums(psf_R), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 20000))


### ### Separate endophytes vs. pathogen
##### Separate symptom vs. no symptom ones###
psf_R_pa <- subset_samples(psf_R, Symptom_condition=="Symptom")
psf_R_en <- subset_samples(psf_R, Symptom_condition=="No_Symptom")
psf_R_bahia<- subset_samples(psf_R, Plant_species=="Paspalum_notatum")
psf_R_bahia_en<- subset_samples(psf_R_bahia, Symptom_condition=="No_Symptom")



## plant
psf_R_m = merge_samples(psf_R, "Plant_species")
sample_data(psf_R_m)$Plant_species <- levels(as.factor(sample_data(psf_R)$Plant_species))

psf_R_m = transform_sample_counts(psf_R_m, function(x) 100 * x/sum(x))

png("barplot_plantall_species.png", width=3000, height=3000,res=350)
title = "fungal class: plant include symptom and non-symptom"
plot_bar(psf_R_m, "Plant_species", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

## barplot endophytes
psf_R_en_m = merge_samples(psf_R_en, "Plant_species")
sample_data(psf_R_en_m)$Plant_species <- levels(as.factor(sample_data(psf_R_en)$Plant_species))

psf_R_en_m = transform_sample_counts(psf_R_en_m, function(x) 100 * x/sum(x))

png("barplot_endo_species.png", width=3000, height=3000,res=350)
title = "endophytes only"
plot_bar(psf_R_en_m, "Plant_species", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()


## barplot pathogen
psf_R_pa_m = merge_samples(psf_R_pa, "Plant_species")
sample_data(psf_R_pa_m)$Plant_species <- levels(as.factor(sample_data(psf_R_pa)$Plant_species))

psf_R_pa_m = transform_sample_counts(psf_R_pa_m, function(x) 100 * x/sum(x))

png("barplot_pathogen_species.png", width=3000, height=3000,res=350)
title = "Pathogen only"
plot_bar(psf_R_pa_m, "Plant_species", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()


## symptom
psf_R_m = merge_samples(psf_R, "Symptom_condition")
#psf_R_mm=merge_taxa(psf_R_m,"Class")
sample_data(psf_R_m)$Symptom_condition <- levels(as.factor(sample_data(psf_R)$Symptom_condition))

psf_R_mm = transform_sample_counts(psf_R_m, function(x) 100 * x/sum(x))

png("symptom_readpercent.png", width=3000, height=3000,res=350)
title = "Symptom_condition"
plot_bar(psf_R_m, "Symptom_condition", fill = "Class", title = title)  + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()


##Location
psf_R_m = merge_samples(psf_R, "Location")
sample_data(psf_R_m)$Location <- levels(as.factor(sample_data(psf_R)$Location))

psf_R_m = transform_sample_counts(psf_R_m, function(x) 100 * x/sum(x))

png("barplot_location.png", width=3000, height=3000,res=350)
title = "Location"
plot_bar(psf_R_m, "Location", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

## neighbor
psf_R_m = merge_samples(psf_R, "Neighbor")
#psf_R_mm=merge_taxa(psf_R_m,"Class")
sample_data(psf_R_m)$Neighbor <- levels(as.factor(sample_data(psf_R)$Neighbor))

psf_R_mm = transform_sample_counts(psf_R_m, function(x) 100 * x/sum(x))

png("barplot_neighbor.png", width=3000, height=3000,res=350)
title = "Neighbor"
plot_bar(psf_R_m, "Neighbor", fill = "Class", title = title)  + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

# Endophytes location
psf_R_en_ml = merge_samples(psf_R_en, "Location")
sample_data(psf_R_en_ml)$Location <- levels(as.factor(sample_data(psf_R_en)$Location))

psf_R_en_ml = transform_sample_counts(psf_R_en_ml, function(x) 100 * x/sum(x))

png("barplot_endophytes_location.png", width=3000, height=3000,res=350)
title = "endophytes only"
plot_bar(psf_R_en_ml, "Location", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

## barplot pathogen
psf_R_pa_m = merge_samples(psf_R_pa, "Location")
sample_data(psf_R_pa_m)$Location <- levels(as.factor(sample_data(psf_R_pa)$Location))

psf_R_pa_m = transform_sample_counts(psf_R_pa_m, function(x) 100 * x/sum(x))

png("barplot_pathogen_location.png", width=3000, height=3000,res=350)
title = "Pathogen only"
plot_bar(psf_R_pa_m, "Location", fill = "Class", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()



######
#### Exam symptom in different plants (NOT WORKING!!!!!)
sample_data(psf_R)$symplant <- paste0(sample_data(psf_R)$Symptom_condition,sample_data(psf_R)$Plant)
psf_symplant = merge_samples(psf_R, "symplant")


# repair factors
sample_data(psf_symplant)$Plant <- levels(as.factor((sample_data(psf_R)$Plant)))[get_variable(psf_symplant,"Plant")]
sample_data(psf_symplant)$Symptom_condition <- levels(as.factor((sample_data(psf_R)$Symptom_condition)))[get_variable(psf_symplant, "Symptom_condition")]
# transform to proportions
psf_symplantP = transform_sample_counts(psf_symplant, function(x) 100 * x/sum(x))

#restroomRgsm19 = prune_taxa(top19otus, restroomRgsmp)
title = "symptom across plant"
p = plot_bar(psf_symplantP, "Symptom_condition",fill = "Class", title = title) + coord_flip() + 
  labs(colour = "Class")
p + facet_wrap(~Plant)
#################


#### PCoA all###
png("Pcoa_all.png", width=3000, height=3000,res=350)
title = "PCoA bray"
psf_ord = ordinate(psf_R, "PCoA", "bray")
p = plot_ordination(psf_R, psf_ord, color = "Plant_species", shape = "Symptom_condition")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
dev.off()

png("Pcoa_all_wrap.png", width=3000, height=3000,res=350)
p+facet_wrap(~Location)
dev.off()


png("Pcoa_pathogen.png", width=3000, height=3000,res=350)
title = "PCoA pathogen"
psf_ord = ordinate(psf_R_pa, "PCoA", "bray")
p = plot_ordination(psf_R_pa, psf_ord, color = "Plant_species",shape="Location")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
dev.off()



png("Pcoa_endophytes.png", width=3000, height=3000,res=350)
title = "PCoA endophytes"
psf_ord = ordinate(psf_R_en, "PCoA", "bray")
p = plot_ordination(psf_R_en, psf_ord, color = "Plant_species")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
dev.off()

png("Pcoa_endophytes_wrap.png", width=3000, height=3000,res=350)
p+facet_wrap(~Location)
dev.off()

### Neighbor Not very informative
title = "PCoA neighbor"
psf_ord = ordinate(psf_R, "PCoA", "bray")
p = plot_ordination(psf_R, psf_ord, color = "Plant",shape="Neighbor")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p

p+facet_wrap(~Neighbor)


### Comparison###
plant_group = get_variable(psf_R, "Genus")
plant_ano = anosim(distance(psf_R, "bray"), plant_group)
plant_ano$signif
plant_ano$statistic

plant_group = get_variable(psf_R_en, "Plant")
plant_ano = anosim(distance(psf_R_en, "bray"), plant_group)
plant_ano$signif
plant_ano$statistic

plant_group = get_variable(psf_R_en, "Genus")
plant_ano = anosim(distance(psf_R_en, "bray"), plant_group)
plant_ano$signif
plant_ano$statistic


loc_group = get_variable(psf_R, "Location")
loc_ano = anosim(distance(psf_R, "bray"), loc_group)
loc_ano$signif
loc_ano$statistic

## Adonis Permutation test
df = as(sample_data(psf_R_en), "data.frame")
d = distance(psf_R_en, "bray")
psf_Radonis = adonis(d ~ Plant + Location= , df)
psf_Radonis

plot(psf_Radonis$aov.tab)


df = as(sample_data(psf_R_en), "data.frame")
d = distance(psf_R_en, "bray")
psf_Radonis = adonis(d ~ Genus , df)
psf_Radonis

plot(psf_Radonis$aov.tab)

df = as(sample_data(psf_R_bahia_en), "data.frame")
d = distance(psf_R_bahia_en, "bray")
psf_Radonis = adonis(d ~ Location , df)
psf_Radonis

## Richness ###

plot_richness(psf_R, x = "Symptom_condition") + geom_boxplot()
plot_richness(psf_R_pa, x = "Plant") + geom_boxplot()
plot_richness(psf_R_pa, x = "Location") + geom_boxplot()

### Networks ##

png("network_Plant.png", width=2000, height=2500,res=250)
ig = make_network(psf_R, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R, color = "Plant_species", shape = "Symptom_condition", line_weight = 0.4, 
             label = NULL)
dev.off()

png("network_Plant_endophytes.png", width=2000, height=2500,res=250)
ig = make_network(psf_R_en, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R_en, color = "Plant_species", line_weight = 0.4, 
             label = NULL)
dev.off()


png("network_Plant_endophytes_bahiagrass.png", width=2000, height=2500,res=250)
ig = make_network(psf_R_en, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R_en, color = "Location", line_weight = 0.4, 
             label = NULL)
dev.off()


png("network_Plant_endophytes_bahiagrass.png", width=2000, height=2500,res=250)
ig = make_network(psf_R_bahia_en, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R_bahia_en, color = "Location", line_weight = 0.4, 
             label = NULL)
dev.off()


###
png("network_neighbor.png", width=2000, height=2500,res=250)
ig = make_network(psf_R, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R, color = "Neighbor", shape="Plant_species", line_weight = 0.4, 
             label = NULL)
dev.off()



ig = make_network(psf_R_pa, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, psf_R_pa, color = "Location", line_weight = 0.4, 
             label = NULL)


mt_psfR_allsym<-mt(psf_R, "Symptom_condition", test="f")
mt_psfR_endo_plant<-mt(psf_R_en,"Plant",test="f")
mt_psfR_pa_plant<-mt(psf_R_pa,"Plant",test="f")

write.table(mt_psfR_allsym,"mt_psfR_allsym_sym.txt", quote=F)
write.table(mt_psfR_endo_plant,"mt_psfR_endo_plant.txt", quote=F)
write.table(mt_psfR_pa_plant,"mt_psfR_pa_plant.txt", quote=F)



## MUST TRY subset_ord_plot ##
#GP <- subset_taxa(GP, Phylum=="Ascomycetes")
#gpca <- ordinate(GP, "PCA")
#p1 = plot_ordination(GP, gpca, "species", color="Class")
#p1





### Rank curves ##



## Various Neighborhood effect add w symptom##

## Nestedness ##


#### Heatmap ##

png("heatmap.png", width=2000, height=2500,res=250)
plot_heatmap(psf_R, "PCoA", "bray","Symptom_condition", sample.order = "Symptom_condition",taxa.label = "Class")
dev.off()


