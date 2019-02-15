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
metaM <-read.delim("/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/metadata_ITS_update.txt")


## Load table for mycotoxin## ZEAR_4_SULF, Ergine and Avenacein_Y_Q1 are only presense/absense
toxin_value<-data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/value.txt"))
#toxin_meta<- data.frame(read.delim("/Users/kc178/Documents/Florida_projects/florencia/mycotoxin_reformat/metadata.txt"))

toxin_value<- apply(toxin_value, 2, function(x){gsub("BLOQ",0.5,x)}) ## Replace cells with peak on HPLC but below detection limit to 0.5
toxin_value[is.na(toxin_value)] <- 0
toxin_value[toxin_value==""]<-0
toxin_value<- data.frame(MycoID=toxin_value[,1], apply(toxin_value[,-1], 2, as.numeric))

toxin_quality<-toxin_value
toxin_quality[toxin_quality>0]<-1
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


psM <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metaM_comV), 
               tax_table(tax))
psM
