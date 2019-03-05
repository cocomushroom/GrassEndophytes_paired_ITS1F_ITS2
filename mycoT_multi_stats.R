data <- sample_data(psMfilter)
otu_table(psMfilter)

data_tv <- metaM_comV[metaM_comV[["MycoID"]]!="No", !names(metaM_comV) %in% c("ITS_ID", "LSU_code")]
data_tq <- metaM_comQ[metaM_comQ[["MycoID"]]!="No",]

data_tv <- unique(data_tv)
unique(data_tv[2:3, ])
dim(data_tv)

t(lapply(unique(data_tv[["MycoID"]]), function(x){
  data_tv[data_tv[["MycoID"]]==x, tox_list]
}))



names(data_tv)
str(data_tv)


table(data_tv[,c("Plant", "Plant_Genus")])
table(data_tv[,c("Plant", "Plant_sci")])
table(data_tv[,c("Plant_Genus", "Plant_species")])
table(data_tv[,c(16,17)])

table(data_tv[,c(16:15)])
table(data_tv[,c(14,16)])
table(data_tv[,c(15,17)])

table(data_tv[,c("Location", "Neighbor")])
table(data_tv[,c("Symptom_condition", "Neighbor")])
table(data_tv[,c("Plant_Genus", "Symptom_condition")])
table(data_tv[,c("Plant_Genus", "Neighbor")])

table(data_tv[,c("Location", "Longitude")])
table(data_tv[,c(18,20)])
table(data_tv[,c(19,20)])

table(data_tv[,c(13,21)])

table(data_tv[, "Symptom_condition"])

tox_list <- names(toxin_value)[-1]
c("MycoID", "ZEAR", "AME", "Emodin", "Alternariol", "AcetylDON",
"Sterigmato", "BEAU", "Citrinin", "ZEAR_4_SULF", "Ergine", "Avenacein_Y_Q1",
"ITS_ID", "Plant", "Plant_Genus", "Plant_sci", "Plant_species",
"Location", "Symptom_condition", "Neighbor", "LSU_code", "Longitude",
"Latitude")
meta_list <- c("MycoID", "Plant", "Plant_species", "Location", "Symptom_condition")
data2 <- data_tv[, c()]

library(MASS)
var_list <- c("Plant", "Plant_Genus", "Plant_sci", "Plant_species",
"Location", "Symptom_condition", "Neighbor")
var <- "Symptom_condition"
var <- "Plant"
var <- "Plant_species"
var <- "Location"

result <- list()
for(var in var_list){
f1 <- glm(formula(sprintf("%s ~ 1", var)), data=data_tv, family=binomial(link = "logit"))
f2 <- glm(formula(sprintf("%s ~ (%s)^2", var, paste(tox_list, collapse = " + "))), data=data_tv, family=binomial(link = "logit"))

summary(f2)
result[[var]] <- stepAIC(f1, scope=list(f1, f2), trace=F)
}




result <- prcomp(data_tv[,tox_list], center=T, scale=T)
summary(result)
par(mfrow=c(1,2))
biplot(result)
plot(result$x[,1:2], col=factor(data_tv[["Symptom_condition"]]))
plot(result$x[,1:2], col=factor(data_tv[["Plant"]]))
plot(result$x[,1:2], col=factor(data_tv[["Location"]]))
names(result)



otu <- otu_table(psMfilter)

result_otu <- list()
for(var in var_list){
  f1 <- glm(data_tv[[var]] ~ 1, family=binomial(link = "logit"))
  f2 <- glm(data_tv[[var]] ~ otu, family=binomial(link = "logit"))
  result_otu[[var]] <- stepAIC(f1, scope=list(f1, f2), trace=F)
}
