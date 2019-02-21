###############################################
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
########### Starting phyloseq##################
setwd("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/phylofig")
library("ape")
library("vegan")
library("tidyr")
library("metacoder")
library(phyloseq)
library("ggplot2")

## Load taxa produced by dada2##
tax30<-readRDS(file = "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/tax30_final.rds")
tax_R<-gsub(".__", "",tax30$tax) ## If having problem with tax rank, skip this
## Load seqtab.nochim produced by dada2##
seqtab.nochim<-readRDS(file = "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/seqtabnochim.rds")

## Not a QIIME2 format
metaM <-read.delim("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/metadata_ITS_update.txt")


## Load table for mycotoxin## ZEAR_4_SULF, Ergine and Avenacein_Y_Q1 are only presense/absense
toxin_value<-data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/value.txt"))
#toxin_meta<- data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/metadata.txt"))

toxin_value<- apply(toxin_value, 2, function(x){gsub("BLOQ",0.5,x)}) ## Replace cells with peak on HPLC but below detection limit to 0.5
toxin_value[is.na(toxin_value)] <- 0
toxin_value[toxin_value==""]<-0
toxin_value<- data.frame(MycoID=toxin_value[,1], apply(toxin_value[,-1], 2, as.numeric))

toxin_quality<-toxin_value
toxin_quality[toxin_quality>0]<-"Yes"
toxin_quality[toxin_quality==0]<-"No"
toxin_quality$MycoID<-toxin_value$MycoID

## Merge with Mycotoxin#
overlap<-toxin_value$MycoID %in% metaM$MycoID
toxin_value_f<-toxin_value[overlap,]

overlap<-toxin_quality$MycoID %in% metaM$MycoID
toxin_quality_f<-toxin_quality[overlap,]

metaM_comV<-merge(toxin_value_f, metaM, by="MycoID",all=TRUE)
row.names(metaM_comV)<-metaM_comV$ITS_ID

metaM_comQ<-merge(toxin_quality_f, metaM, by="MycoID",all=TRUE)
row.names(metaM_comQ)<-metaM_comQ$ITS_ID

#merge(metaM, toxin_quality,by="MycoID")

## note this still have non fungi, toxin quality first
psM <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metaM_comQ), 
               tax_table(tax_R))

psMfilter <-  prune_samples(sample_names(psM)[!sample_data(psM)$MycoID =="No"],psM) #sample_data(psM)$MycoID =="No" # also worked?

psMfilterF <- subset_taxa(psMfilter, Kingdom=="Fungi")
nreads = sort(sample_sums(psMfilter)) # if need the lowest read number for a sample to determine the rarification number
nreads


#### Create Rarified table ### 
set.seed(1111)
psMfilterF_R = rarefy_even_depth(psMfilterF, sample.size = nreads[[1]])
saveRDS(psMfilterF_R, "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psMfilterF_R.rds")


#### Metacoder###
psMfilterF_R<-readRDS(file="/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/psMfilterF_R.rds") # already rarified in Phyloseq
psMfilterF_Rp <- parse_phyloseq(psMfilterF_R)
print(psMfilterF_Rp) #pay attention to different tables and names

psMfilterF_Rp$data$otu_table <- zero_low_counts(psMfilterF_Rp, data = "otu_table", min_count = 2)

no_reads <- rowSums(psMfilterF_Rp$data$otu_table[, psMfilterF_Rp$data$sample_data$sample_id]) == 0
sum(no_reads)

psMfilterF_Rp <- filter_obs(psMfilterF_Rp, data = "tax_data", ! no_reads, drop_taxa = TRUE)
print(psMfilterF_Rp)

psMfilterF_Rp$data$tax_abund <- calc_taxon_abund(psMfilterF_Rp, "otu_table",
                                          cols = psMfilterF_Rp$data$sample_data$sample_id)

psMfilterF_Rp$data$tax_occ <- calc_n_samples(psMfilterF_Rp, "tax_abund", groups = psMfilterF_Rp$data$sample_data$Citrinin)

print(psMfilterF_Rp)


heat_tree(psMfilterF_Rp,
          node_size = rowMeans(cbind(Yes, No)),
          node_color =  Yes - No,
          node_label = taxon_names,
          tree_label = taxon_names,
          node_color_range = diverging_palette(),
          #node_color_range = c("black", "#FFFFFF"), #use for change main color
          #node_color_interval = c(-0.1, 0.1),
          #edge_color_interval = c(-0.1, 0.1),
          #edge_color_range = c("black", "#FFFFFF"),
          #node_label_max = 200,
          layout = "da", initial_layout = "re",
          title = "Citrinin",
          node_size_axis_label = "Average proportion",
          output_file = "metacoder_all_toxinQ_Citrinin.pdf")


psMfilterF_Rp %>%
  taxa::filter_taxa(taxon_names == "Ascomycota", subtaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  taxa::filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  heat_tree(
            node_size = rowMeans(cbind(Yes, No)),
            node_color =  Yes - No,
            node_label = taxon_names,
            tree_label = taxon_names,
            node_color_range = diverging_palette(),
            #node_color_range = c("black", "#FFFFFF"), #use for change main color
            #node_color_interval = c(-0.1, 0.1),
            #edge_color_interval = c(-0.1, 0.1),
            #edge_color_range = c("black", "#FFFFFF"),
            #node_label_max = 200,
            layout = "da", initial_layout = "re",
            title = "ZEAR",
            node_size_axis_label = "Average proportion",
            output_file = "metacoder_all_toxinQ_ZEAR_Asco_Genus.pdf")

png("zear4_asco_genus.png", width=9000, height=8000,res=350)
psMfilterF_Rp %>%
  taxa::filter_taxa(taxon_names == "Ascomycota", subtaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  taxa::filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% # subset to the class rank # Need to change if fix taxonomy in the future
  heat_tree(
    node_size = rowMeans(cbind(Yes, No)),
    node_color =  Yes - No,
    node_label = taxon_names,
    tree_label = taxon_names,
    node_color_range = diverging_palette(),
    #node_color_range = c("black", "#FFFFFF"), #use for change main color
    #node_color_interval = c(-0.1, 0.1),
    #edge_color_interval = c(-0.1, 0.1),
    #edge_color_range = c("black", "#FFFFFF"),
    #node_label_max = 200,
    layout = "da", initial_layout = "re",
    title = "Zear4",
    node_size_axis_label = "Average proportion") #,output_file = "metacoder_all_toxinQ_ZEAR4_Basi_Genus.pdf"
dev.off()