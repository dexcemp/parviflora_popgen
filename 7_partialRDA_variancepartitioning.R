##Parviflora group

#Extract dbMEMs
library(adespatial)
parvifloracoords <- read.csv("parviflora_ss_coords.csv", row.names = 1)
parviflora_dbmem <- dbmem(parvifloracoords)
parviflora_dbmem <- as.data.frame(parviflora_dbmem)

#Run RDA with dbMEMs
parvicomplex.rda.dbmem <- rda(snp_cleaned ~., data = parviflora_dbmem)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_cleaned ~ 1, data = parviflora_dbmem),
                      scope = formula(parvicomplex.rda.dbmem),
                      direction = "forward",
                      R2scope = TRUE,
                      pstep = 9999,
                      trace = FALSE)

# Check the new model with forward-selected variables
fwd.sel$call

#New model
parvicomplex.rda.dbmem.signif <- rda(snp_cleaned ~ MEM2 + MEM4 + MEM3 + MEM6 + MEM1 + 
                                       MEM16 + MEM7 + MEM14 + MEM9 + MEM5 + MEM10 + MEM8 + MEM11, data = parviflora_dbmem) #Top 3: MEM2, MEM4, and MEM3

#Significance testing
anova.cca(parvicomplex.rda.dbmem.signif, step = 1000, by = "term") #MEM11 not significant

#New model; significant dbems only
parvicomplex.rda.dbmem.signif <- rda(snp_cleaned ~ MEM2 + MEM4 + MEM3 + MEM6 + MEM1 + 
                                       MEM16 + MEM7 + MEM14 + MEM9 + MEM5 + MEM10 + MEM8, data = parviflora_dbmem)

#Merge environmental data and dbMEMs
parviflora_ss_combined <- cbind(env.z, parviflora_dbmem)

#Run partial RDA (full model conditioned on all significant dbMEMs)
parvicomplex.rda.spacecontrolled <- rda(snp_cleaned ~ BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_15 +      
                                          BIOCLIM_8 + BIOCLIM_9 + GTOPO30_ASPECT_reduced + 
                                          GTOPO30_ELEVATION + GTOPO30_SLOPE_reduced + ISRICSOILGRIDS_new_average_claypercent_reduced +      
                                          ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced +      
                                          ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + 
                                          Condition(MEM2 + MEM4 + MEM3 + MEM6 + MEM1 + MEM16 + MEM7 + MEM14 + MEM9 + MEM5 + MEM10 + MEM8), data = parviflora_ss_combined)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(parvicomplex.rda.spacecontrolled) 

#Significance testing
anova.cca(parvicomplex.rda.spacecontrolled, step = 1000) #still significant
anova.cca(parvicomplex.rda.spacecontrolled, step = 1000, by = "term") #by term
summary(parvicomplex.rda.spacecontrolled)

#Variance partitioning
#Load data
envi <-read.csv("parviflora_ss_envipoints.csv", header=T, row.names = 1)
envi <- as.data.frame(envi)

#Run PCA
envi_PCA <- prcomp(envi, scale. = TRUE)

#Scores (observations)
enviscores <- as.data.frame(envi_PCA$x)

#Run RDA with PCs
rdaPC <- rda(snp_cleaned ~., data = enviscores)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_cleaned ~ 1, data = enviscores), 
                      scope = formula(rdaPC), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#New model
rdaPC.signif <- rda(snp_cleaned ~ PC3 + PC1 + PC2 + PC7 + PC5 + PC12 + 
                      PC9 + PC8 + PC18 + PC16 + PC27 + PC26 + PC4 + PC24 + PC19 + 
                      PC10 + PC20 + PC6 + PC25 + PC14 + PC11 + PC21 + PC17 + PC13 + 
                      PC28 + PC22 + PC15, data = enviscores) 

#Significance testing
anova.cca(rdaPC.signif, step = 1000, by = "term") #Top 3: 3, 1, 2

#Merge environmental PCs and dbMEMs
combined <- cbind(enviscores, parviflora_dbmem)

#Run variance partitioning
#Environment only
rdaPC.envi <- rda(snp_cleaned ~ PC3 + PC1 + PC2 + Condition (MEM2 + MEM4 + MEM3), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.envi)

#Significance testing
anova.cca(rdaPC.envi, step = 1000)
summary(rdaPC.envi)

#Space only
rdaPC.space <- rda(snp_cleaned ~ MEM2 + MEM4 + MEM3 + Condition (PC3 + PC1 + PC2), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.space)

#Significance testing
anova.cca(rdaPC.space, step = 1000)
summary(rdaPC.space)

#Both environment and space
rdaPC.both <- rda(snp_cleaned ~ PC3 + PC1 + PC2 + MEM2 + MEM4 + MEM3, data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.both)

#Significance testing
anova.cca(rdaPC.both, step = 1000)
summary(rdaPC.both)


##H. puberula
#Extract dbMEMs
library(adespatial)
puberulacoords <- read.csv("puberula_coords.csv", row.names = 1)
puberula_dbmem <- dbmem(puberulacoords)
puberula_dbmem <- as.data.frame(puberula_dbmem)

#Run RDA with dbMEMs
puberula.rda.dbmem <- rda(snppub_cleaned ~., data = puberula_dbmem)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snppub_cleaned ~ 1, data = puberula_dbmem),
                      scope = formula(puberula.rda.dbmem), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#New model
puberula.rda.dbmem.signif <- rda(snppub_cleaned ~ MEM1 + MEM2 + MEM3 + MEM4, data = puberula_dbmem) #Top 3: MEM1, MEM2, and MEM3

#Significance testing
anova.cca(puberula.rda.dbmem.signif, step = 1000, by = "term") #all significant

#Merge environmental data and dbMEMs
puberula_combined <- cbind(env_pub.scaled, puberula_dbmem)

##Run partial RDA (full model conditioned on all significant dbMEMs)
puberula.rda.spacecontrolled <- rda(snppub_cleaned ~ BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_15 +
                                      BIOCLIM_2 + GTOPO30_ASPECT_reduced + GTOPO30_SLOPE_reduced +
                                      ISRICSOILGRIDS_new_average_bulkdensity_reduced +
                                      ISRICSOILGRIDS_new_average_phx10percent_reduced +
                                      ISRICSOILGRIDS_new_average_sandpercent_reduced + Condition(MEM1 + MEM2 + MEM3 + MEM4), data = puberula_combined)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(puberula.rda.spacecontrolled)

#Significance testing
anova.cca(puberula.rda.spacecontrolled, step = 1000) #still significant
anova.cca(puberula.rda.spacecontrolled, step = 1000, by = "term") #by term
summary(puberula.rda.spacecontrolled)

#Variance partitioning
#Load data
envi <-read.csv("puberula_envipoints.csv", header=T, row.names = 1)
envi <- as.data.frame(envi)

#Run PCA
envi_PCA <- prcomp(envi, scale. = TRUE)

#Scores (observations)
enviscores <- as.data.frame(envi_PCA$x)

#Run RDA with PCs
rdaPC <- rda(snppub_cleaned ~., data = enviscores)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(snppub_cleaned ~ 1, data = enviscores),
                      scope = formula(rdaPC),
                      direction = "forward",
                      R2scope = TRUE,
                      pstep = 9999,
                      trace = FALSE)

#Check the new model with forward-selected variables
fwd.sel$call

#New model
rdaPC.signif <- rda(snppub_cleaned ~ PC28 + PC19 + PC22 + PC25 + PC7 + 
                      PC17 + PC16 + PC6 + PC12 + PC10 + PC23 + PC1, data = enviscores) 

#Significance testing
anova.cca(rdaPC.signif, step = 1000, by = "term") #Top 3: PC28, PC19, PC22


#Merge environmental PCs and dbMEMs
combined <- cbind(enviscores, puberula_dbmem)

#Run variance partitioning

#Environment only
rdaPC.envi <- rda(snppub_cleaned ~ PC28 + PC19 + PC22 + Condition (MEM1 + MEM2 + MEM3), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.envi)

#Significance testing
anova.cca(rdaPC.envi, step = 1000)
summary(rdaPC.envi)

#Space only
rdaPC.space <- rda(snppub_cleaned ~ MEM1 + MEM2 + MEM3 + Condition (PC28 + PC19 + PC22), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.space)

#Significance testing
anova.cca(rdaPC.space, step = 1000)
summary(rdaPC.space)

#Both environment and space
rdaPC.both <- rda(snppub_cleaned ~ PC28 + PC19 + PC22 + MEM1 + MEM2 + MEM3, data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.both)

#Significance testing
anova.cca(rdaPC.both, step = 1000)
summary(rdaPC.both)


##H. missouriensis
#Extract dbMEMs
library(adespatial)
missocoords <- read.csv("misso_coords.csv", row.names = 1)
misso_dbmem <- dbmem(missocoords)
misso_dbmem <- as.data.frame(misso_dbmem)

#Run RDA with dbMEMs
misso.rda.dbmem <- rda(snpmisso_cleaned ~., data = misso_dbmem)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpmisso_cleaned ~ 1, data = misso_dbmem), 
                      scope = formula(misso.rda.dbmem),
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE)

#Check the new model with forward-selected variables
fwd.sel$call #removed MEM5

#New model
misso.rda.dbmem.signif <- rda(snpmisso_cleaned ~ MEM2 + MEM4 + MEM3 + MEM1 + 
                                MEM6, data = misso_dbmem) #Top 3: MEM2, MEM4, and MEM3

anova.cca(misso.rda.dbmem.signif, step = 1000, by = "term") #all significant

#Merge environmental data and dbMEMs
misso_combined <- cbind(env_misso.scaled, misso_dbmem)

#Run partial RDA (full model conditioned on all significant dbMEMs)
misso.rda.spacecontrolled <- rda(snpmisso_cleaned ~ BIOCLIM_1 + BIOCLIM_10 + BIOCLIM_12 +
                                   BIOCLIM_15 + BIOCLIM_8 + GTOPO30_ASPECT_reduced + GTOPO30_ELEVATION +
                                   ISRICSOILGRIDS_new_average_claypercent_reduced +
                                   ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced +
                                   ISRICSOILGRIDS_new_average_phx10percent_reduced + Condition(MEM2 + MEM4 + MEM3 + MEM1 + MEM6), data = misso_combined)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(misso.rda.spacecontrolled)

#Significance testing
anova.cca(misso.rda.spacecontrolled, step = 1000) #still significant
anova.cca(misso.rda.spacecontrolled, step = 1000, by = "term") #by term
summary(misso.rda.spacecontrolled)

#Variance partitioning
#Load data
envi <-read.csv("misso_envipoints.csv", header=T, row.names = 1)
envi <- as.data.frame(envi)

#Run PCA
envi_PCA <- prcomp(envi, scale. = TRUE)

# Scores (observations)
enviscores <- as.data.frame(envi_PCA$x)

#Run RDA with PCs
rdaPC <- rda(snpmisso_cleaned ~., data = enviscores)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpmisso_cleaned ~ 1, data = enviscores), 
                      scope = formula(rdaPC), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#New model
rdaPC.signif <- rda(snpmisso_cleaned ~ PC23 + PC5 + PC10 + PC7 + PC1 + 
                      PC3 + PC6 + PC8 + PC9 + PC24 + PC11 + PC4 + PC22 + PC28 + 
                      PC16 + PC17 + PC26 + PC21 + PC12, data = enviscores) 

#Significance testing
anova.cca(rdaPC.signif, step = 1000, by = "term") #Top 3: PC23, PC5, PC10

#Merge environmental PCs and dbMEMs
combined <- cbind(enviscores, misso_dbmem)

#Run variance partitioning
#Environment only
rdaPC.envi <- rda(snpmisso_cleaned ~ PC23 + PC5 + PC10 + Condition (MEM2 + MEM4 + MEM3), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.envi) 

#Significance testing
anova.cca(rdaPC.envi, step = 1000)
summary(rdaPC.envi)

#Space only
rdaPC.space <- rda(snpmisso_cleaned ~ MEM2 + MEM4 + MEM3 + Condition (PC23 + PC5 + PC10), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.space) 

#Significance testing
anova.cca(rdaPC.space, step = 1000)
summary(rdaPC.space)

#Both environment and space
rdaPC.both <- rda(snpmisso_cleaned ~ PC23 + PC5 + PC10 + MEM2 + MEM4 + MEM3, data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.both) 

#Significance testing
anova.cca(rdaPC.both, step = 1000)
summary(rdaPC.both)


##H. parviflora var. parviflora
#Extract dbMEMs
library(adespatial)
varparvicoords <- read.csv("varparvi_coords.csv", row.names = 1)
varparvi_dbmem <- dbmem(varparvicoords)
varparvi_dbmem <- as.data.frame(varparvi_dbmem)

#Run RDA with dbMEMs
varparvi.rda.dbmem <- rda(snp_varparvi_cleaned ~., data = varparvi_dbmem)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_varparvi_cleaned ~ 1, data = varparvi_dbmem), 
                      scope = formula(varparvi.rda.dbmem), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call #removed MEM4

#New model
varparvi.rda.dbmem.signif <- rda(snp_varparvi_cleaned ~ MEM3 + MEM2 + MEM1 + MEM6 + 
                                   MEM7 + MEM5, data = varparvi_dbmem) #Top 3:MEM3, MEM2, and MEM1

#Significance testing
anova.cca(varparvi.rda.dbmem.signif, step = 1000, by = "term") #all significant

#Merge environmental data and dbMEMs
varparvi_combined <- cbind(env_varparvi.scaled, varparvi_dbmem)

#Run partial RDA (full model conditioned on all significant dbMEMs)
varparvi.rda.spacecontrolled <- rda(snp_varparvi_cleaned ~ BIOCLIM_1 + BIOCLIM_10 + BIOCLIM_12 +
                                      BIOCLIM_15 + BIOCLIM_2 + BIOCLIM_8 + BIOCLIM_9 + GTOPO30_ASPECT_reduced +
                                      GTOPO30_SLOPE_reduced + ISRICSOILGRIDS_new_average_claypercent_reduced +
                                      ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced +
                                      ISRICSOILGRIDS_new_average_phx10percent_reduced +
                                      ISRICSOILGRIDS_new_average_sandpercent_reduced + Condition(MEM3 + MEM2 + MEM1 + MEM6 + MEM7 + MEM5), data = varparvi_combined)

# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(varparvi.rda.spacecontrolled) 

#Significance testing
anova.cca(varparvi.rda.spacecontrolled, step = 1000) 
anova.cca(varparvi.rda.spacecontrolled, step = 1000, by = "term") #by term
summary(varparvi.rda.spacecontrolled)

#Variance partitioning
#Load data
envi <-read.csv("varparvi_envipoints.csv", header=T, row.names = 1)
envi <- as.data.frame(envi)

#Run PCA
envi_PCA <- prcomp(envi, scale. = TRUE)

# Scores (observations)
enviscores <- as.data.frame(envi_PCA$x)

#Run RDA with PCs
rdaPC <- rda(snp_varparvi_cleaned ~., data = enviscores)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_varparvi_cleaned ~ 1, data = enviscores),
                      scope = formula(rdaPC), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 
#Check the new model with forward-selected variables
fwd.sel$call

#New model
rdaPC.signif <- rda(snp_varparvi_cleaned ~ PC9 + PC25 + PC20 + PC7 + 
                      PC2 + PC27 + PC13 + PC8 + PC5 + PC15 + PC3 + PC22 + PC16 + 
                      PC12 + PC14, data = enviscores) 

anova.cca(rdaPC.signif, step = 1000, by = "term") #Top 3: PC9, PC25, PC20


#Merge environemntal PCs and dbMEMs
combined <- cbind(enviscores, varparvi_dbmem)

#Run variance partitioning
#Environment only
rdaPC.envi <- rda(snp_varparvi_cleaned ~ PC9 + PC25 + PC20 + Condition (MEM3 + MEM2 + MEM1), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.envi)

#Significance testing
anova.cca(rdaPC.envi, step = 1000)
summary(rdaPC.envi)

#Space only
rdaPC.space <- rda(snp_varparvi_cleaned ~ MEM3 + MEM2 + MEM1 + Condition(PC9 + PC25 + PC20), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.space)

#Significance testing
anova.cca(rdaPC.space, step = 1000)
summary(rdaPC.space)

#Both environment and space
rdaPC.both <- rda(snp_varparvi_cleaned ~ PC9 + PC25 + PC20 + MEM3 + MEM2 + MEM1, data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.both)

#Significance testing
anova.cca(rdaPC.both, step = 1000)
summary(rdaPC.both)


##H. parviflora var. saurensis
#Extract dbMEMs
library(adespatial)
varsaurcoords <- read.csv("varsaur_coords.csv", row.names = 1)
varsaur_dbmem <- dbmem(varsaurcoords)
varsaur_dbmem <- as.data.frame(varsaur_dbmem)

#Run RDA with dbMEMs
varsaur.rda.dbmem <- rda(snpsaur_cleaned ~., data = varsaur_dbmem)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpsaur_cleaned ~ 1, data = varsaur_dbmem), 
                      scope = formula(varsaur.rda.dbmem), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call #same, just MEM1

#New model
varsaur.rda.dbmem.signif <- rda(snpsaur_cleaned ~ MEM1, data = varsaur_dbmem) #just MEM1

#Significance testing
anova.cca(varsaur.rda.dbmem.signif, step = 1000, by = "term")

#Merge environmental data and dbMEMs
varsaur_combined <- cbind(env_saur.scaled, varsaur_dbmem)

#Run partial RDA (full model conditioned on all significant dbMEMs)
varsaur.rda.spacecontrolled <- rda(snpsaur_cleaned ~ BIOCLIM_15 + GTOPO30_SLOPE_reduced +
                                     ISRICSOILGRIDS_new_average_phx10percent_reduced + Condition (MEM1), data = varsaur_combined)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(varsaur.rda.spacecontrolled) 

#Significance testing
anova.cca(varsaur.rda.spacecontrolled, step = 1000) #still significant
anova.cca(varsaur.rda.spacecontrolled, step = 1000, by = "term") #by term
summary(varsaur.rda.spacecontrolled)

#Variance partitioning
#Load data
envi <-read.csv("varsaur_envipoints.csv", header=T, row.names = 1)
envi <- as.data.frame(envi)

#Run PCA
envi_PCA <- prcomp(envi, scale. = TRUE)

# Scores (observations)
enviscores <- as.data.frame(envi_PCA$x)

#Run RDA with PCs
rdaPC <- rda(snpsaur_cleaned ~., data = enviscores)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpsaur_cleaned ~ 1, data = enviscores), 
                      scope = formula(rdaPC), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#New model
rdaPC.signif <- rda(snpsaur_cleaned ~ PC2, data = enviscores) 

#Significance testing
anova.cca(rdaPC.signif, step = 1000, by = "term") #top1: PC2


#Merge environmental PCs and dbMEMs
combined <- cbind(enviscores, varsaur_dbmem)

#Run variance partitioning

#Environment only
rdaPC.envi <- rda(snpsaur_cleaned ~ PC2 + Condition (MEM1), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.envi) 

#Significance testing
anova.cca(rdaPC.envi, step = 1000)
summary(rdaPC.envi)

#Space only
rdaPC.space <- rda(snpsaur_cleaned ~ MEM1 + Condition(PC2), data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.space) 

#Significance testing
anova.cca(rdaPC.space, step = 1000)
summary(rdaPC.space)

#Both environment and space
rdaPC.both <- rda(snpsaur_cleaned ~ PC2 + MEM1, data = combined) 

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rdaPC.both)

#Significance testing
anova.cca(rdaPC.both, step = 1000)
summary(rdaPC.both)
