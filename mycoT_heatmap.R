library(gplots)
library(RColorBrewer)
data_tv <- metaM_comV[metaM_comV[["MycoID"]]!="No", !names(metaM_comV) %in% c("ITS_ID", "LSU_code")]

data_tv <- unique(data_tv)
dim(data_tv)
length(unique(data_tv[["MycoID"]]))

tox_list <- names(toxin_value)[-1]

data_m <- as.matrix(data_tv[,tox_list])


# distMethods <- c(Euclidean="euclidean", Maximum='maximum', Manhattan='manhattan', Canberra='canberra', Binary='binary', Minkowski='minkowski')
#
# hclustMethods <- c(Complete="complete", Single="single", Average="average", Mcquitty="mcquitty", Median="median", Centroid="centroid", Ward.D="ward.D", Ward.D2= "ward.D2")

cPal255 <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n=255)

row.names(data_m) <- data_tv[["Plant_Genus"]]

for(var in c("Plant_Genus", "Location", "Neighbor")){
  png(sprintf("mycotoxin_%s.png", var), width=800, height=800)
  row.names(data_m) <- data_tv[[var]]
  data_log <- log(data_m+1)
  heatmap.2(data_log,
    col=cPal255,
    distfun=function(x) {dist(x, method="euclidean")},
    hclustfun=function(x) {hclust(x, method="complete")},
    dendrogram="both", Rowv=TRUE, Colv=TRUE,
    symm=F, symkey=F, symbreaks=F,
    key=T,  trace="none", density.info="none",
    lwid=c(1, 4), lhei=c(1, 4),
    margin=c(10,10))
  dev.off()
}



data_norm <- apply(data_m, 2, function(x){x/max(x)})
heatmap.2(data_norm,
  col=cPal255,
  distfun=function(x) {dist(x, method="euclidean")},
  hclustfun=function(x) {hclust(x, method="single")},
  dendrogram="both", Rowv=TRUE, Colv=TRUE,
  symm=F, symkey=F, symbreaks=F,
  key=T, trace="none", density.info="none",
  lwid=c(1, 4), lhei=c(1, 4),
  margin=c(10,10))
