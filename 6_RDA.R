#Extract climate data
library(raster)

#Make list of file paths for each layer
clim.list <- dir("/mnt/Heuheu/climate_data/Saxifragales_all_layers_30s", full.names=T, pattern = "\\.tif$")  

#Stack the layers into a single object
clim.layer <-  stack(clim.list)  

#Crop the climate data layers to just the area of interest. I defined an extent with minimum and maximum longitude and minimum and maximum latitude.
extent <- c(-93, -80, 33, 39)
clim.layer.crop <- crop(clim.layer, extent)

#Load our sample coordinates and extract climate data at each of the points.
sample.coord <-read.table("samplecoords_parviflora_ss.txt", header=T, stringsAsFactors=F)
sample.coord

#Define the spatial projection system that the points are in (usually WGS84)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  
sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('longitude','latitude')], proj4string=CRS(crs.wgs), data=sample.coord)

#Extract the data for each point (projection of climate layer and coordinates must match)
clim.points <- extract(clim.layer.crop, sample.coord.sp)  

#Combines the sample coordinates with the climate data points
clim.points <- cbind(sample.coord, clim.points)  

#Save the table for downstream analyses
write.table(clim.points, "parviflora_ss_clim.points", sep="\t", quote=F, row.names=F) 

#Calculate Pearson correlation
#Read csv containing env values
combined <- read.csv("parviflora_ss_corrtest.csv", header = TRUE)

#Read columns containing numbers as a matrix
df <- combined[ -1 ]
combmatrix <- as.matrix(df)

#Calculate correlations
cors <- cor(combmatrix, method="pearson", use="everything")

##Get p-value significance
pvals <- matrix(NA, ncol(combmatrix), ncol(combmatrix))

for(i in 1:ncol(combmatrix)) {
  for(j in 1:ncol(combmatrix)) {
    test <- cor.test(combmatrix[, i], combmatrix[, j])
    pvals[i, j] <- test$p.value
  }
}

#Plot heat map
heatmap(cors, margins=c(10,10), revC=T, col=rev(heat.colors(256)))

#Write results; numbers above 0.8 represent high correlation between the variables
write.csv(cors, "parviflora_ss.pearson.csv")
write.csv(pvals, "parviflora_ss_pearson_pvals.csv")


#Redundancy Analyses (RDA)

##Parviflora group
#Read snp and env file
snp <- read.table("parviflora_ss_snp.forR", header = T, row.names = 1)
env <- read.csv("parviflora_ss_uncorrelated.csv", row.names = 1)

#Check for NAs
summary(snp)
anyNA(snp)
any(is.nan(as.matrix(snp)))
any(is.infinite(as.matrix(snp)))

#Mean imputation
snp_cleaned <- apply(snp, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

#Standardizing the environmental variables
library(vegan)
#Scale and center variables
env.z <- decostand(env, method = "standardize")

#Variables are now centered around a mean of 0
round(apply(env.z, 2, mean), 1)

#Scaled to have a standard deviation of 1
apply(env.z, 2, sd)

#Run the RDA
parvicomplex.rda <- rda(snp_cleaned ~., data = env.z)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_cleaned ~ 1, data = env.z), 
                      scope = formula(parvicomplex.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#Write our new model
parvicomplex.rda.signif <- rda(snp_cleaned ~ ISRICSOILGRIDS_new_average_sandpercent_reduced + 
                                 BIOCLIM_8 + BIOCLIM_15 + ISRICSOILGRIDS_new_average_claypercent_reduced + 
                                 ISRICSOILGRIDS_new_average_phx10percent_reduced + GTOPO30_ELEVATION + 
                                 ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + 
                                 GTOPO30_SLOPE_reduced + BIOCLIM_9 + BIOCLIM_12 + BIOCLIM_1 + 
                                 GTOPO30_ASPECT_reduced, data = env.z)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(parvicomplex.rda.signif)

#Significance testing
anova.cca(parvicomplex.rda.signif, step = 1000)
anova.cca(parvicomplex.rda.signif, step = 1000, by = "term") #by term
summary(parvicomplex.rda.signif)

#Visualization

#Extract RDA scores
site_scores <- vegan::scores(parvicomplex.rda.signif,
                             display="sites",
                             choices=c(1,2),
                             scaling=1)

env_scores <- vegan::scores(parvicomplex.rda.signif,
                            display="bp",
                            choices=c(1,2),
                            scaling=1)

#Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$sample <- rownames(site_df)
env_df <- as.data.frame(env_scores)
env_df$variable <- rownames(env_df)

#Add population metadata
metadata <- read.table("parviflora_ss_popdata.txt", header=TRUE)
metadata <- as.data.frame(metadata)
site_df$Pop <- metadata$Pop

#Add proportion explained by each RDA axis
perc <- round(100*(summary(parvicomplex.rda.signif)$cont$importance[2, 1:2]), 2)

#Scale arrows so they appear clearly
env_df$RDA1 <- env_df$RDA1 * 2
env_df$RDA2 <- env_df$RDA2 * 2

#Plot in ggplot2
library(ggplot2)
ggplot(site_df, aes(x=RDA1, y=RDA2, fill=Pop)) +
  geom_point(shape = 21, size=3, color = "black") +
  scale_fill_manual(values = c(
    "puberula" = "#00946a",
    "missouriensis" = "#facc14",
    "parviflora" = "#305cde",
    "saurensis" = "#ec5800")) +
  geom_segment(data=env_df,
               aes(x=0, y=0, xend=RDA1, yend=RDA2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.2,"cm")),
               color="black") +
  geom_text(data=env_df,
            aes(x=RDA1, y=RDA2, label=variable),
            size = 2,
            inherit.aes = FALSE,
            color="black") +
  xlab(paste0("RDA1 (", perc[1], "%)"))+ 
  ylab(paste0("RDA2 (", perc[2], "%)"))+ 
  theme_bw()

##H. puberula
#Read snp and env file
snp_pub <- read.table("puberula_snp.forR", header = T, row.names = 1)
env_pub <- read.csv("puberula_uncorrelated.csv", row.names = 1)

#Check NAs
anyNA(snp_pub)
any(is.nan(as.matrix(snp_pub)))
any(is.infinite(as.matrix(snp_pub)))

#Mean imputation
snppub_cleaned <- apply(snp_pub, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
anyNA(snppub_cleaned)

#Standardizing the environmental variables
library(vegan)
#Scale and center variables
env_pub.scaled <- decostand(env_pub, method = "standardize")

#Variables are now centered around a mean of 0
round(apply(env_pub.scaled, 2, mean), 1)

#Scaled to have a standard deviation of 1
apply(env_pub.scaled, 2, sd)

#Run the RDA
puberula.rda <- rda(snppub_cleaned ~., data = env_pub.scaled)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snppub_cleaned ~ 1, data = env_pub.scaled), 
                      scope = formula(puberula.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE)

#Check the new model with forward-selected variables
fwd.sel$call

#Write our new model
puberula.rda.signif <- rda(snppub_cleaned ~ BIOCLIM_1 + BIOCLIM_2 + BIOCLIM_12 + 
                             BIOCLIM_15 + GTOPO30_SLOPE_reduced + GTOPO30_ASPECT_reduced + 
                             ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_bulkdensity_reduced + 
                             ISRICSOILGRIDS_new_average_phx10percent_reduced, data = env_pub.scaled)
#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(puberula.rda.signif)

#Significance testing
anova.cca(puberula.rda.signif, step = 1000)
anova.cca(puberula.rda.signif, step = 1000, by = "term") #by term
summary(puberula.rda.signif)

#Visualization
#Extract RDA scores
site_scores <- vegan::scores(puberula.rda.signif,
                             display="sites",
                             choices=c(1,2),
                             scaling=1)

env_scores <- vegan::scores(puberula.rda.signif,
                            display="bp",
                            choices=c(1,2),
                            scaling=1)

#Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$sample <- rownames(site_df)
env_df <- as.data.frame(env_scores)
env_df$variable <- rownames(env_df)

#Add population metadata
metadata <- read.table("puberula_popdata.txt", header=TRUE)
metadata <- as.data.frame(metadata)
site_df$Pop <- metadata$Pop

##Add proportion explained by each RDA axis
perc <- round(100*(summary(puberula.rda.signif)$cont$importance[2, 1:2]), 2)

#Scale arrows so they appear clearly
env_df$RDA1 <- env_df$RDA1 * 2
env_df$RDA2 <- env_df$RDA2 * 2

#Plot in ggplot2
library(ggplot2)
ggplot(site_df, aes(x=RDA1, y=RDA2, fill=Pop)) +
  geom_point(shape = 21, size=4, color = "black") +
  scale_fill_manual(values = c(
    "North" = "#00946a",
    "South" = "#8bcf00")) +
  geom_segment(data=env_df,
               aes(x=0, y=0, xend=RDA1, yend=RDA2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.2,"cm")),
               color="black") +
  geom_text(data=env_df,
            aes(x=RDA1, y=RDA2, label=variable),
            size = 2,
            inherit.aes = FALSE,
            color="black") +
  xlab(paste0("RDA1 (", perc[1], "%)"))+ 
  ylab(paste0("RDA2 (", perc[2], "%)"))+ 
  theme_bw()


##H. missouriensis
#Read snp and env file
snp_misso <- read.table("misso_snp.forR", header = T, row.names = 1)
env_misso <- read.csv("misso_uncorrelated.csv", row.names = 1)

#Check for NAs
anyNA(snp_misso)
any(is.nan(as.matrix(snp_misso)))
any(is.infinite(as.matrix(snp_misso)))

#Mean imputation
snpmisso_cleaned <- apply(snp_misso, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
anyNA(snpmisso_cleaned)

#Standardizing the environmental variables
library(vegan)
#Scale and center variables
env_misso.scaled <- decostand(env_misso, method = "standardize")

#Variables are now centered around a mean of 0
round(apply(env_misso.scaled, 2, mean), 1)

#Scaled to have a standard deviation of 1
apply(env_misso.scaled, 2, sd)

#Run the RDA
misso.rda <- rda(snpmisso_cleaned ~., data = env_misso.scaled)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpmisso_cleaned ~ 1, data = env_misso.scaled), 
                      scope = formula(misso.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE)

#Check the new model with forward-selected variables
fwd.sel$call

#Write our new model
misso.rda.signif <- rda(snpmisso_cleaned ~ BIOCLIM_15 + BIOCLIM_12 + ISRICSOILGRIDS_new_average_phx10percent_reduced + 
                          BIOCLIM_1 + BIOCLIM_10 + GTOPO30_ELEVATION + BIOCLIM_8 + 
                          ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + 
                          ISRICSOILGRIDS_new_average_claypercent_reduced + GTOPO30_ASPECT_reduced, 
                        data = env_misso.scaled)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(misso.rda.signif)

#Significance testing
anova.cca(misso.rda.signif, step = 1000)
anova.cca(misso.rda.signif, step = 1000, by = "term") #by term
summary(misso.rda.signif)

#Visualization
#Extract RDA scores
site_scores <- vegan::scores(misso.rda.signif,
                             display="sites",
                             choices=c(1,2),
                             scaling=1)

env_scores <- vegan::scores(misso.rda.signif,
                            display="bp",
                            choices=c(1,2),
                            scaling=1)

#Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$sample <- rownames(site_df)
env_df <- as.data.frame(env_scores)
env_df$variable <- rownames(env_df)

#Add population metadata
metadata <- read.table("misso_popdata.txt", header=TRUE)
metadata <- as.data.frame(metadata)
site_df$Pop <- metadata$Pop

#Add proportion explained by each RDA axis
perc <- round(100*(summary(misso.rda.signif)$cont$importance[2, 1:2]), 2)

#Scale arrows so they appear clearly
env_df$RDA1 <- env_df$RDA1 * 2
env_df$RDA2 <- env_df$RDA2 * 2

#Plot in ggplot2
library(ggplot2)
ggplot(site_df, aes(x=RDA1, y=RDA2, fill=Pop)) +
  geom_point(shape = 21, size=4, color = "black") +
  scale_fill_manual(values = c(
    "a" = "#facc14",
    "b" = "#fffb00",
    "c" = "#a98600",
    "d" = "#e9d700",
    "e"= "#8c8645")) +
  geom_segment(data=env_df,
               aes(x=0, y=0, xend=RDA1, yend=RDA2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.2,"cm")),
               color="black") +
  geom_text(data=env_df,
            aes(x=RDA1, y=RDA2, label=variable),
            size = 2,
            inherit.aes = FALSE,
            color="black") +
  xlab(paste0("RDA1 (", perc[1], "%)"))+ 
  ylab(paste0("RDA2 (", perc[2], "%)"))+ 
  theme_bw()


##H. parviflora var. parviflora
#Read snp and env file
snp_varparvi <- read.table("varparvi_snp.forR", header = T, row.names = 1)
env_varparvi <- read.csv("varparvi_uncorrelated.csv", row.names = 1)

#Check for NAs
anyNA(snp_varparvi)
any(is.nan(as.matrix(snp_varparvi)))
any(is.infinite(as.matrix(snp_varparvi)))

#Mean imputation
snp_varparvi_cleaned <- apply(snp_varparvi, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
anyNA(snp_varparvi_cleaned)

#Standardizing the environmental variables
library(vegan)
#Scale and center variables
env_varparvi.scaled <- decostand(env_varparvi, method = "standardize")

#Variables are now centered around a mean of 0
round(apply(env_varparvi.scaled, 2, mean), 1)

#Scaled to have a standard deviation of 1
apply(env_varparvi.scaled, 2, sd)

#Run the RDA
varparvi.rda <- rda(snp_varparvi_cleaned ~., data = env_varparvi.scaled)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snp_varparvi_cleaned ~ 1, data = env_varparvi.scaled),
                      scope = formula(varparvi.rda ), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Check the new model with forward-selected variables
fwd.sel$call

#Write our new model
varparvi.rda.signif <- rda(snp_varparvi_cleaned ~ BIOCLIM_15 + ISRICSOILGRIDS_new_average_phx10percent_reduced + 
                             BIOCLIM_9 + BIOCLIM_10 + BIOCLIM_1 + BIOCLIM_8 + ISRICSOILGRIDS_new_average_sandpercent_reduced + 
                             BIOCLIM_12 + GTOPO30_ASPECT_reduced + ISRICSOILGRIDS_new_average_claypercent_reduced + 
                             ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + 
                             GTOPO30_SLOPE_reduced, data = env_varparvi.scaled)

#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(varparvi.rda.signif)

#Significance testing
anova.cca(varparvi.rda.signif, step = 1000)
anova.cca(varparvi.rda.signif, step = 1000, by = "term") #by term
summary(varparvi.rda.signif)

#Visuualization
#Extract RDA scores
site_scores <- vegan::scores(varparvi.rda,
                             display="sites",
                             choices=c(1,2),
                             scaling=1)

env_scores <- vegan::scores(varparvi.rda,
                            display="bp",
                            choices=c(1,2),
                            scaling=1)

#Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$sample <- rownames(site_df)
env_df <- as.data.frame(env_scores)
env_df$variable <- rownames(env_df)

#Add population metadata
metadata <- read.table("varparvi_popdata.txt", header=TRUE)
metadata <- as.data.frame(metadata)
site_df$Pop <- metadata$Pop

#Add proportion explained by each RDA axis
perc <- round(100*(summary(varparvi.rda)$cont$importance[2, 1:2]), 2)

#Scale arrows so they appear clearly
env_df$RDA1 <- env_df$RDA1 * 2
env_df$RDA2 <- env_df$RDA2 * 2

#Plot in ggplot2
library(ggplot2)
ggplot(site_df, aes(x=RDA1, y=RDA2, fill=Pop)) +
  geom_point(shape = 21, size=4, color = "black") +
  scale_fill_manual(values = c(
    "N" = "#305cde",
    "S" = "#6fd0ee")) +
  geom_segment(data=env_df,
               aes(x=0, y=0, xend=RDA1, yend=RDA2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.2,"cm")),
               color="black") +
  geom_text(data=env_df,
            aes(x=RDA1, y=RDA2, label=variable),
            size = 2,
            inherit.aes = FALSE,
            color="black") +
  xlab(paste0("RDA1 (", perc[1], "%)"))+ 
  ylab(paste0("RDA2 (", perc[2], "%)"))+ 
  theme_bw()


##H. parviflora var. saurensis
#Read snp and env file
snp_saur <- read.table("varsaur_snp.forR", header = T, row.names = 1)
env_saur <- read.csv("varsaur_uncorrelatednew.csv", row.names = 1)

#Check for NAs
anyNA(snp_saur)
any(is.nan(as.matrix(snp_saur)))
any(is.infinite(as.matrix(snp_saur)))

#Mean imputation
snpsaur_cleaned <- apply(snp_saur, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
anyNA(snpsaur_cleaned)

#Standardizing the environmental variables
library(vegan)
#Scale and center variables
env_saur.scaled <- decostand(env_saur, method = "standardize")

#Variables are now centered around a mean of 0
round(apply(env_saur.scaled, 2, mean), 1)

#Scaled to have a standard deviation of 1
apply(env_saur.scaled, 2, sd)

#Run the RDA
varsaur.rda <- rda(snpsaur_cleaned ~., data = env_saur.scaled)

#Forward selection of variables:
fwd.sel <- ordiR2step(rda(snpsaur_cleaned ~ 1, data = env_saur.scaled), 
                      scope = formula(varsaur.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 9999,
                      trace = FALSE) 

#Run new model
varsaur.rda.signif <- rda(snpsaur_cleaned ~ BIOCLIM_15 + ISRICSOILGRIDS_new_average_phx10percent_reduced, 
                          data = env_saur.scaled)


#Check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(varsaur.rda.signif)

#Significance testing
anova.cca(varsaur.rda.signif, step = 1000)
anova.cca(varsaur.rda.signif, step = 1000, by = "term") #by term
summary(varsaur.rda.signif)

#Visualization
#Extract RDA scores
site_scores <- vegan::scores(varsaur.rda.signif,
                             display="sites",
                             choices=c(1,2),
                             scaling=1)

env_scores <- vegan::scores(varsaur.rda.signif,
                            display="bp",
                            choices=c(1,2),
                            scaling=1)

#Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$sample <- rownames(site_df)
env_df <- as.data.frame(env_scores)
env_df$variable <- rownames(env_df)

#Add population metadata
metadata <- read.table("varsaur_popdata.txt", header=TRUE)
metadata <- as.data.frame(metadata)
site_df$Pop <- metadata$Pop

#Add proportion explained by each RDA axis
perc <- round(100*(summary(varsaur.rda.signif)$cont$importance[2, 1:2]), 2)

#Scale arrows so they appear clearly
env_df$RDA1 <- env_df$RDA1 * 2
env_df$RDA2 <- env_df$RDA2 * 2

#Plot in ggplot2
library(ggplot2)
ggplot(site_df, aes(x=RDA1, y=RDA2, fill=Pop)) +
  geom_point(shape = 21, size=4, color = "black") +
  scale_fill_manual(values = c(
    "a" = "#ec5800",
    "b" = "#ffa000",
    "c" = "#b46202")) +
  geom_segment(data=env_df,
               aes(x=0, y=0, xend=RDA1, yend=RDA2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.2,"cm")),
               color="black") +
  geom_text(data=env_df,
            aes(x=RDA1, y=RDA2, label=variable),
            size = 2,
            inherit.aes = FALSE,
            color="black") +
  xlab(paste0("RDA1 (", perc[1], "%)"))+ 
  ylab(paste0("RDA2 (", perc[2], "%)"))+ 
  theme_bw()
