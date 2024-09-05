#=======================================#
#                                       #
# Global Model                          #
#                                       #
# Linguistics cooperation               #
# Bernd Lenzner                         #
# bernd.lenzner@univie.ac.at            #
# 2023                                  #
#                                       #
#=======================================#
datum <- Sys.Date()


#===========================================#
# Load relevant packages ----
## Package names of the ones we need
packages <- c("stringr","ncdf4","tictoc","tidyverse", "sf","data.table","GGally","lme4","glmmTMB","cowplot","sjPlot","tinytable", "corrplot", "ggcorrplot", "ggpubr", "DHARMa", "gamlss", "performance") 

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Load packages
invisible(lapply(packages, library, character.only = TRUE))



#===========================================#
# Directory and paths ----
drive <- "E"
path.datasets <- paste0(drive, ':/Datasets') # path for directory with predictor datasets
path.proj <- paste0(drive, ':/Publication_Projects/2022_Linguistics_Biodiv_Language_Threat') # path for directory with code
path.data <- paste0(path.proj, '/Data') # path for directory with analysis datasets
path.fig <- paste0(path.proj, '/Figures/2024-09-Final-figures') # path for directory with analysis datasets
path.code <- paste0(path.proj, '/R-Scripts') # path for directory with analysis datasets

# Functions ----
source(paste0(path.code,"/F01-plotting-predictor-correlations.R"))


# Shapefile ----
## Final region selection (countries GADM level 0) ----
gadm.data.path <- paste0(path.datasets, "/Regions/GADM/GADM0-simpl-0.1.shp")
gadm0 <- st_read(gadm.data.path) # national simplified
shape <- gadm0 %>%
  st_wrap_dateline() %>% # avoid banding
  st_transform(crs = "+proj=moll")

sf::sf_use_s2(FALSE) # disable s2 proccessing to run st_union() function
shape.proj <- shape

# Full datasets ----
load(paste0(path.data, "/2024-04-22_Full-Analysis-Dataset.Rdata"))
load(paste0(path.data, "/2024-07-09-RedList-cat-ALL-by-region.RData")) # RedList categories by group and region

# ISO-code
ISO <- read_delim(paste0(path.datasets, "/Regions/ISO/ISO_inc_Regions.csv"))
ISO.sub <- ISO %>%
  select(`alpha-3`, region, `sub-region`)

#=======================================================================================#
# Build datasets for each response variable ----
## Drop correlated predictors ----
data.global <- data.full %>%
  select(!c(dem.range, aridity.sd, crop.area.sum, pop.2015, struct_div, pop.dens.2015,crop.dens.2019, hdi.pc.2015.mean, traveltime.mean1, traveltime.mean9,
            Alpha.2.code))

#=======================================================================================#
## Data subsets and and variable transformation ----
### LDI and RLI animals combined ----
df.biocultural <- data.global %>%
  select(!c(RLI.no.amph:RLI.rept, total_languages, avg_score_languages)) %>%
  mutate(Ind.mean = rowMeans(.[,c("LDI", "RLI.animal")])) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(gdp.pc.2015.mean),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)
  )  %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))


### Language Diversity Index ----
df.linguistic <- data.global %>%
  select(!c(RLI.animal:RLI.rept, total_languages, avg_score_languages)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)
  )  %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))


### Red List Index - Animals ----
df.biodiversity <- data.global %>%
  select(!c(RLI.no.amph:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)
         ) %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3")) 


### Red List Index - Amphibians ----
df.rli.amph <- data.global %>%
  select(!c(RLI.animal:RLI.no.rept, RLI.birds:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)) %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))


### Red List Index - Birds ----
df.rli.bird <- data.global %>%
  select(!c(RLI.animal:RLI.amph, RLI.mamm:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)) %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))

### Red List Index - mammals ----
df.rli.mamm <- data.global %>%
  select(!c(RLI.animal:RLI.birds, RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)) %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))


### Red List Index - reptiles ----
df.rli.rept <- data.global %>%
  select(!c(RLI.animal:RLI.mamm, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)),
         dem.s = scale(dem.sd),
         aridity.s = scale(aridity.mean),
         temp.s = scale(temp.mean),
         gdp.s = scale(log(gdp.pc.2015.mean + 0.1)),
         schooling.s = scale(`age.2015_15+`),
         traveltime.large.s = scale(log(traveltime.median1 + 0.1)),
         traveltime.small.s = scale(log(traveltime.median9 + 0.1)),
         urban.s = scale(log(urban + 0.1)),
         intense.s = scale(log(`intense-prod` + 0.1)),
         mixed.s =  scale(log(`mixed-cult` + 0.1)),
         occ.time.full.s = scale(occ.time)) %>% 
  left_join(ISO.sub, by = c("Alpha.3.code" = "alpha-3"))


### Red List Index - Animals (no amphibians) ----
df.rli.no.amph <- data.global %>%
  select(!c(RLI.animal, RLI.no.bird:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)))


### Red List Index - Animals (no birds) ----
df.rli.no.bird <- data.global %>%
  select(!c(RLI.animal:RLI.no.amph, RLI.no.mamm:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)))


### Red List Index - Animals (no mammals) ----
df.rli.no.mamm <- data.global %>%
  select(!c(RLI.animal:RLI.no.bird, RLI.no.rept:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)))


### Red List Index - Animals (no reptiles) ----
df.rli.no.rept <- data.global %>%
  select(!c(RLI.animal:RLI.no.mamm, RLI.amph:RLI.rept, LDI, total_languages, avg_score_languages, struct_dist, genetic_div)) %>%
  drop_na(.) %>%
  mutate(area.s = scale(log(area + 0.1)),
         samp.bias.s = scale(log(samp.bias.mean + 0.1)))



## Save correlation of indicator subsets ----
ggsave(filename = paste0(path.fig, "/", datum, "-Figure-CORRELATION-indicator-subsets-Language.png"),
       cowplot::plot_grid(plot.predictor.correlation("lang.RLI"),
                          plot.predictor.correlation("LDI"),
                          plot.predictor.correlation("RLI.animal"),
                          ncol = 1),
       width = 7, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-Figure-CORRELATION-indicator-subsets-animals.pdf"),
       cowplot::plot_grid(plot.predictor.correlation("RLI.amph"),
                          plot.predictor.correlation("RLI.bird"),
                          plot.predictor.correlation("RLI.mamm"),
                          plot.predictor.correlation("RLI.rept"), ncol = 2),
       width = 10, height = 14, unit = "in", bg = "white")


ggsave(filename = paste0(path.fig, "/", datum, "-Figure-CORRELATION-indicator-subsets-animals1.png"),
       cowplot::plot_grid(plot.predictor.correlation("RLI.amph"),
                          plot.predictor.correlation("RLI.bird"), ncol = 2),
       width = 10, height = 7, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-Figure-CORRELATION-indicator-subsets-animals2.png"),
       cowplot::plot_grid(plot.predictor.correlation("RLI.mamm"),
                          plot.predictor.correlation("RLI.rept"), ncol = 2),
       width = 10, height = 7, unit = "in", bg = "white")




#==============================================================================#
#==============================================================================#
#==============================================================================#
# Hotspot analysis ----
## Account for area and sampling bias from data ----
### LDI/RLI -indicator (mean) ----
m.hotspot.biocultural <- gamlss(Ind.mean ~ area.s, 
                           data = df.biocultural,
                           family = BE(),
                           method = mixed())

summary(m.hotspot.biocultural)

### Language ----
df.linguistic$LDI2 <- scales::rescale(df.linguistic$LDI, from=c(0, 1), to = c(0.0001, 0.9999)) # rescale variable to exclude values of 1

m.hotspot.linguistic <- gamlss(LDI2~ area.s, 
                      data = df.linguistic,
                      family = BE(),
                      method = mixed())

summary(m.hotspot.linguistic)

### Animals ----
m.hotspot.biodiversty <- gamlss(RLI.animal ~ area.s * samp.bias.s, 
                        data = df.biodiversity,
                        family = BE(),
                        method = mixed())

summary(m.hotspot.biodiversty)

### Amphibians ----
df.rli.amph$RLI.amph2 <- scales::rescale(df.rli.amph$RLI.amph, from=c(0, 1), to = c(0.0001, 0.9999)) # rescale variable to exclude values of 1
m.hotspot.amph <- gamlss(RLI.amph2 ~ area.s * samp.bias.s, 
                      data = df.rli.amph,
                      family = BE(),
                      method = mixed())

### Birds ----
m.hotspot.bird <- gamlss(RLI.birds ~ area.s * samp.bias.s, 
                      data = df.rli.bird,
                      family = BE(),
                      method = mixed())

### Mammals ----
m.hotspot.mamm <- gamlss(RLI.mamm ~ area.s * samp.bias.s, 
                      data = df.rli.mamm,
                      family = BE(),
                      method = mixed())

### Reptiles ----
df.rli.rept$RLI.rept2 <- scales::rescale(df.rli.rept$RLI.rept, from=c(0, 1), to = c(0.0001, 0.9999)) # rescale variable to exclude values of 1
m.hotspot.rept <- gamlss(RLI.rept2 ~ area.s * samp.bias.s, 
                      data = df.rli.rept,
                      family = BE(),
                      method = mixed())

### No amphibians ----
m.hotspot.no.amph <- gamlss(RLI.no.amph ~ area.s * samp.bias.s, 
                         data = df.rli.no.amph,
                         family = BE(),
                         method = mixed())
### No birds ----
m.hotspot.no.bird <- gamlss(RLI.no.bird ~ area.s * samp.bias.s, 
                         data = df.rli.no.bird,
                         family = BE(),
                         method = mixed())

### No mammals ----
m.hotspot.no.mamm <- gamlss(RLI.no.mamm ~ area.s * samp.bias.s, 
                         data = df.rli.no.mamm,
                         family = BE(),
                         method = mixed())

### No reptiles ----
m.hotspot.no.rept <- gamlss(RLI.no.rept ~ area.s * samp.bias.s, 
                         data = df.rli.no.rept,
                         family = BE(),
                         method = mixed())

## Save hotspot models ----
save(m.hotspot.biocultural, m.hotspot.linguistic, m.hotspot.biodiversty, m.hotspot.amph, m.hotspot.bird, m.hotspot.mamm, m.hotspot.rept, file = paste0(path.data,"/",datum, "-Hotspot-Models.RData"))


#### Model Deviance - beta model ----
data.frame(Model = c("LDI/RLI [mean]","Languages (LDI)", "Animals (RLI)", "Amphibians (RLI)", "Birds (RLI)", "Mammals (RLI)", "Reptiles (RLI)"),
           R.squared = c(round(Rsq(m.hotspot.biocultural),2),
                        round(Rsq(m.hotspot.linguistic),2),
                        round(Rsq(m.hotspot.biodiversty),2),
                        round(Rsq(m.hotspot.amph),2),
                        round(Rsq(m.hotspot.bird),2),
                        round(Rsq(m.hotspot.mamm),2),
                        round(Rsq(m.hotspot.rept),2)))


## Check for autocorrelation ----
### ACF Plots ----
#### Check for autocorrelation - beta model
acf(m.hotspot.biocultural$residuals, type="correlation") # LDI/RLI Indikator (mean)
acf(m.hotspot.linguistic$residuals, type="correlation") # LDI
acf(m.hotspot.biodiversty$residuals, type="correlation") # RLI (animals)
acf(m.hotspot.amph$residuals, type="correlation") # RLI (amphibians)
acf(m.hotspot.bird$residuals, type="correlation") # RLI (birds)
acf(m.hotspot.mamm$residuals, type="correlation") # RLI (mammals)
acf(m.hotspot.rept$residuals, type="correlation") # RLI (reptiles)

### Durbin-Watson Test ----
### Beta Model: DW between 1.5 & 2.5 implies no autocorrelation
check_autocorrelation(m.hotspot.biocultural)
check_autocorrelation(m.hotspot.linguistic)
check_autocorrelation(m.hotspot.biodiversty)
check_autocorrelation(m.hotspot.amph)
check_autocorrelation(m.hotspot.bird)
check_autocorrelation(m.hotspot.mamm)
check_autocorrelation(m.hotspot.rept)

#=======================================================================================#
#=======================================================================================#
#=======================================================================================#
# Residual analysis ----
## Add residuals to dataset ----
### Beta regression residuals ----
df.biocultural$Ind.m.resid.beta <- residuals(m.hotspot.biocultural)
df.linguistic$LDI.resid.beta <- residuals(m.hotspot.linguistic)
df.biodiversity$RLI.animal.resid.beta <- residuals(m.hotspot.biodiversty)
df.rli.amph$RLI.amph.resid.beta <- residuals(m.hotspot.amph)
df.rli.bird$RLI.bird.resid.beta <- residuals(m.hotspot.bird)
df.rli.mamm$RLI.mamm.resid.beta <- residuals(m.hotspot.mamm)
df.rli.rept$RLI.rept.resid.beta <- residuals(m.hotspot.rept)

df.rli.no.amph$RLI.no.amph.resid.beta <- residuals(m.hotspot.no.amph)
df.rli.no.bird$RLI.no.bird.resid.beta <- residuals(m.hotspot.no.bird)
df.rli.no.mamm$RLI.no.mamm.resid.beta <- residuals(m.hotspot.no.mamm)
df.rli.no.rept$RLI.no.rept.resid.beta <- residuals(m.hotspot.no.rept)

## Save dataframes ----
save(df.biocultural, df.linguistic, df.biodiversity, df.rli.amph, df.rli.bird, df.rli.mamm, df.rli.rept, file = paste0(path.data,"/",datum, "-Dataframes.RData"))

#=======================================================================================#
## Residual correlation across patterns ----
data.corr.pattern <- data.global %>%
  select(English.short.name:Numeric) %>%
  left_join(df.biocultural %>% select(Alpha.3.code, Ind.m.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.linguistic %>% select(Alpha.3.code, LDI.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.biodiversity %>% select(Alpha.3.code, RLI.animal.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.amph %>% select(Alpha.3.code, RLI.amph.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.bird %>% select(Alpha.3.code, RLI.bird.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.mamm %>% select(Alpha.3.code, RLI.mamm.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.rept %>% select(Alpha.3.code, RLI.rept.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.no.amph %>% select(Alpha.3.code, RLI.no.amph.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.no.bird %>% select(Alpha.3.code, RLI.no.bird.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.no.mamm %>% select(Alpha.3.code, RLI.no.mamm.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.rli.no.rept %>% select(Alpha.3.code, RLI.no.rept.resid.beta), by = "Alpha.3.code")
  

dat.cor.ind <- round(cor(na.omit(data.corr.pattern %>% 
                    select(!c(English.short.name, Alpha.3.code, Numeric))), method = "pearson"),2)

round(cor(na.omit(data.corr.pattern %>% 
                                   select(c(LDI.resid.beta, RLI.animal.resid.beta)))),2)


dat.cor.ind[3,c(4:11)] <- NA
dat.cor.ind[4,c(3,8:11)] <- NA
dat.cor.ind[5,c(3,8:11)] <- NA
dat.cor.ind[6,c(3,8:11)] <- NA
dat.cor.ind[7,c(3,8:11)] <- NA
dat.cor.ind[8,c(3,5:7, 9:11)] <- NA
dat.cor.ind[9,c(3,4,6:7, 8,10:11)] <- NA
dat.cor.ind[10,c(3,4,5,7, 8:9, 11)] <- NA
dat.cor.ind[11,c(3,4:6, 8:10)] <- NA

p.corr.indicators <- ggcorrplot(dat.cor.ind, 
                          #hc.order = TRUE, 
                          type = "lower",
                          title = "Indicator correlations (residuals)",
                          outline.col = "white",
                          ggtheme = ggplot2::theme_bw,
                          colors = c("#6D9EC1", "white", "#E46726"),
                          lab = TRUE,
                          lab_size = 8,
                          tl.cex = 20,
                          insig = "blank") +
  font("title", size = 20, face = "bold") +
  theme(legend.position = "none") +
  font("legend.title", size = 16) +
  font("legend.text", size = 16)

### Save correlation of different indicator residuals ----
ggsave(filename = paste0(path.fig, "/", datum, "-Figure-CORRELATION-indicator-residuals.pdf"),
       p.corr.indicators,
       width = 15, height = 15, unit = "in", bg = "white")






#==============================================================================#
#==============================================================================#
# Driver models - beta regression mixed models ----
#### Biocultural diversity threat ----
m.biocult.beta <- gamlss(Ind.mean ~ area.s + samp.bias.s +
                              dem.s + aridity.s + 
                              traveltime.small.s + urban.s + intense.s + mixed.s + schooling.s + #gdp.s + 
                              occ.time.full.s+random(factor(region)),
                            family = BE(),
                            method = mixed(),
                            data = df.biocultural)
summary(m.biocult.beta)

#### Linguistic diversity threat ----
m.linguistic.beta <- gamlss(LDI2 ~ area.s + 
                              dem.s + aridity.s + 
                              traveltime.small.s + urban.s + intense.s + mixed.s + schooling.s + #gdp.s + 
                              occ.time.full.s+random(factor(region)),
                            family = BE(),
                            method = mixed(),
                            data = df.linguistic)
summary(m.linguistic.beta)

#### Animal diversity threat ----
m.animal.beta <- gamlss(RLI.animal ~ area.s * samp.bias.s + dem.s + aridity.s + 
                                 traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + #schooling.s +
                                 occ.time.full.s+random(factor(region)),
                               family = BE(),
                               method = mixed(),
                               data = df.biodiversity)
summary(m.animal.beta)

#### Amphibian diversity threat ----
m.amph.beta <- gamlss(RLI.amph2 ~ area.s * samp.bias.s + dem.s + aridity.s + 
                               traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + #schooling.s +
                               occ.time.full.s+random(factor(region)),
                             family = BE(),
                             method = mixed(),
                             data = df.rli.amph)
summary(m.amph.beta)

#### Bird diversity threat ----
m.bird.beta <- gamlss(RLI.birds ~ area.s * samp.bias.s + dem.s + aridity.s + 
                               traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + #schooling.s +
                               occ.time.full.s+random(factor(region)),
                             family = BE(),
                             method = mixed(),
                             data = df.rli.bird)
summary(m.bird.beta)

#### Mammal diversity threat ----
m.mamm.beta <- gamlss(RLI.mamm ~ area.s * samp.bias.s + dem.s + aridity.s + 
                               traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + #schooling.s +
                               occ.time.full.s+random(factor(region)),
                             family = BE(),
                             method = mixed(),
                             data = df.rli.mamm)
summary(m.mamm.beta)

#### Reptile diversity threat ----
m.rept.beta <- gamlss(RLI.rept2 ~ area.s * samp.bias.s + dem.s + aridity.s + 
                               traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + #schooling.s +
                               occ.time.full.s+random(factor(region)),
                             family = BE(),
                             method = mixed(),
                             data = df.rli.rept)
summary(m.rept.beta)


## Check model autocorrelation ----
acf(m.biocult.beta$residuals, type="correlation") # Biocultural diversity
acf(m.linguistic.beta$residuals, type="correlation") # Linguistic diversity
acf(m.animal.beta$residuals, type="correlation") # Biodiversity

acf(m.amph.beta$residuals, type="correlation") # Ampibians
acf(m.bird.beta$residuals, type="correlation") # Birds
acf(m.mamm.beta$residuals, type="correlation") # Mammals
acf(m.rept.beta$residuals, type="correlation") # Reptiles

##### Durbin-Watson Test for autocorrelation (performance package)
check_autocorrelation(m.biocult.beta) # Biocultural diversity
check_autocorrelation(m.linguistic.beta) # Linguistic diversity
check_autocorrelation(m.animal.beta) # Biodiversity

check_autocorrelation(m.amph.beta) # Ampibians
check_autocorrelation(m.bird.beta) # Birds
check_autocorrelation(m.mamm.beta) # Mammals
check_autocorrelation(m.rept.beta) # Reptiles


#### Diagnostic plots ----
plot(m.biocult.beta) # Biocultural diversity
plot(m.linguistic.beta) # Linguistic diversity
plot(m.animal.beta) # Biodiversity

plot(m.amph.beta) # Amphibians
plot(m.bird.beta) # Birds
plot(m.mamm.beta) # Mammals
plot(m.rept.beta) # Reptiles


#### Model R-squared ----
data.frame(
  Model = 
    c("Biocultural", "Linguistic", "Animals", "Amphibians", "Birds", "Mammals", "Reptiles"),
  R.squared = 
  round(c(Rsq(m.biocult.beta),
  Rsq(m.linguistic.beta),
  Rsq(m.animal.beta),
  
  Rsq(m.amph.beta),
  Rsq(m.bird.beta),
  Rsq(m.mamm.beta),
  Rsq(m.rept.beta)), 2)
  )


#### Coefficient plot for driver analysis ----
m.drivers <- list()
m.drivers[[1]] <- m.biocult.beta
m.drivers[[2]] <- m.linguistic.beta
m.drivers[[3]] <- m.animal.beta


p.drivers <- plot_models(m.drivers,
                            title = "",
                            show.values = F,
                            auto.label = T,
                            vline.color = "darkgrey",
                            p.shape = T,
                            spacing = 0.6,
                            #colors = c("#0243b3", "#199ebd", "#b21807", "#7d12e4", "#087224", "#675128", "#ee5533"),
                            colors = c("#1485f0", "#f014e6", "#f09914"),
                            m.labels = c("RLI (biocultural)", "RLI (linguistic)", "RLI (biodiversity)"),
                            axis.lim = c(0.3, 1.5),
                            legend.title = "",
                            dot.size	= 6,
                            line.size = 1)  + 
  scale_x_discrete(limits=rev(c("area.s", "samp.bias.s", "area.s:samp.bias.s", "dem.s", "aridity.s", "traveltime.small.s", "urban.s", "intense.s", "mixed.s", "schooling.s", "gdp.s", "occ.time.full.s")),
                   labels = rev(c("Area", "Sampling effort", "Area x Sampling effort","Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "Years of schooling", "GDPpc", "Occupation time"))) + 
  theme_bw() +
  theme(text = element_text(size = 26),
        legend.position = "bottom") + 
  ylim(0.3, 1.5)



ggsave(filename = paste0(path.fig, "/", datum, "-Figure-Drivers-Beta-regression.pdf"),
       p.drivers,
       width = 15, height = 10, unit = "in", bg = "white")


#








#==============================================================================#
# Threat categories by empire ----
RL.cat.all

emp.cats <- RL.cat.all %>%
  select(!c(occ.time.min, occ.time.max, occ.time)) %>% 
  rename_with(~ gsub("occ.time.", "", .x, fixed = TRUE)) %>%
  pivot_longer(!c(GID_0 : n.colonizers), , names_to = "empire") %>%
  mutate(value = case_when(value > 0 ~ empire)) %>%
  select(-empire) %>%
  distinct()

View(emp.cats)

## Calculate threat categories by empire and taxon ----
emp.threats <- emp.cats %>%
  group_by(value, taxon) %>% 
  summarise(across(EN:EX, sum)) %>%
  filter(!is.na(value)) %>%
  pivot_longer(!c(value, taxon), names_to = "threat.cat", values_to = "n.spp")
  

### Grouped barplot
emp.threats$threat.cat <- factor(emp.threats$threat.cat, levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EX"))

ggplot(emp.threats, aes(fill=value, y=n.spp, x=threat.cat)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~taxon) +
  theme_bw() +
  xlab("")

ggplot(emp.threats, aes(x=value, y=n.spp, fill=threat.cat)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~taxon) +
  theme_bw() +
  xlab("")

#==============================================================================#











































#=======================================================================================#
#=======================================================================================#
#=======================================================================================#
#=======================================================================================#
# Driver analysis ----
## Residual Models ----
### LDI/RLI Indicator (mean) ----
m.resid.LDI.RLI.m.beta <- lm(Ind.m.resid.beta ~ dem.s + aridity.s + 
                               traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                               occ.time.full.s,
                             df.biocultural)

summary(m.resid.LDI.RLI.m.beta)

### Language ----
m.resid.lang.beta <- lm(LDI.resid.beta ~ dem.s + aridity.s + 
                          traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                          occ.time.full.s,
                        df.linguistic)

summary(m.resid.lang.beta)

#### Animals ----
m.resid.animal.beta <- lm(RLI.animal.resid.beta ~ dem.s + aridity.s + 
                            traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                            occ.time.full.s,
                          df.biodiversity)

summary(m.resid.animal.beta)

#### Amphibians ----
m.resid.amph.beta <- lm(RLI.amph.resid.beta ~ dem.s + aridity.s + 
                          traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                          occ.time.full.s,
                        df.rli.amph)

summary(m.resid.amph.beta)

#### Birds ----
m.resid.bird.beta <- lm(RLI.bird.resid.beta ~ dem.s + aridity.s + 
                          traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                          occ.time.full.s,
                        df.rli.bird)

summary(m.resid.bird.beta)

#### Mammals ----
m.resid.mamm.beta <- lm(RLI.mamm.resid.beta ~ dem.s + aridity.s + 
                          traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                          occ.time.full.s,
                        df.rli.mamm)

summary(m.resid.mamm.beta)

#### Reptiles ----
m.resid.rept.beta <- lm(RLI.rept.resid.beta ~ dem.s + aridity.s + 
                          traveltime.small.s + urban.s + intense.s + mixed.s + gdp.s + schooling.s +
                          occ.time.full.s,
                        df.rli.rept)

summary(m.resid.rept.beta)

## Model visualizations ----
### LDI/RLI indicator (mean) ----
#### Diagnostic plots ----
p.diag.lang.rli.m1 <- plot_model(m.resid.LDI.RLI.m.beta, type = "diag")[[1]]
p.diag.lang.rli.m2 <- plot_model(m.resid.LDI.RLI.m.beta, type = "diag")[[2]]
p.diag.lang.rli.m3 <- plot_model(m.resid.LDI.RLI.m.beta, type = "diag")[[3]]
p.diag.lang.rli.m4 <- plot_model(m.resid.LDI.RLI.m.beta, type = "diag")[[4]]

### LDI ----
#### Diagnostic plots ----
p.diag.lang1 <- plot_model(m.resid.lang.beta, type = "diag")[[1]]
p.diag.lang2 <- plot_model(m.resid.lang.beta, type = "diag")[[2]]
p.diag.lang3 <- plot_model(m.resid.lang.beta, type = "diag")[[3]]
p.diag.lang4 <- plot_model(m.resid.lang.beta, type = "diag")[[4]]

### RLI - animals ----
#### Diagnostic plots ----
p.diag.animal1 <- plot_model(m.resid.animal.beta, type = "diag")[[1]]
p.diag.animal2 <- plot_model(m.resid.animal.beta, type = "diag")[[2]]
p.diag.animal3 <- plot_model(m.resid.animal.beta, type = "diag")[[3]]
p.diag.animal4 <- plot_model(m.resid.animal.beta, type = "diag")[[4]]

### RLI - amphibians ----
#### Diagnostic plots ----
p.diag.amph1 <- plot_model(m.resid.amph.beta, type = "diag")[[1]]
p.diag.amph2 <- plot_model(m.resid.amph.beta, type = "diag")[[2]]
p.diag.amph3 <- plot_model(m.resid.amph.beta, type = "diag")[[3]]
p.diag.amph4 <- plot_model(m.resid.amph.beta, type = "diag")[[4]]

### RLI - Birds ----
#### Diagnostic plots ----
p.diag.bird1 <- plot_model(m.resid.bird.beta, type = "diag")[[1]]
p.diag.bird2 <- plot_model(m.resid.bird.beta, type = "diag")[[2]]
p.diag.bird3 <- plot_model(m.resid.bird.beta, type = "diag")[[3]]
p.diag.bird4 <- plot_model(m.resid.bird.beta, type = "diag")[[4]]

### RLI - mammals ----
#### Diagnostic plots ----
p.diag.mamm1 <- plot_model(m.resid.mamm.beta, type = "diag")[[1]]
p.diag.mamm2 <- plot_model(m.resid.mamm.beta, type = "diag")[[2]]
p.diag.mamm3 <- plot_model(m.resid.mamm.beta, type = "diag")[[3]]
p.diag.mamm4 <- plot_model(m.resid.mamm.beta, type = "diag")[[4]]

### RLI - reptiles ----
#### Diagnostic plots ----
p.diag.rept1 <- plot_model(m.resid.rept.beta, type = "diag")[[1]]
p.diag.rept2 <- plot_model(m.resid.rept.beta, type = "diag")[[2]]
p.diag.rept3 <- plot_model(m.resid.rept.beta, type = "diag")[[3]]
p.diag.rept4 <- plot_model(m.resid.rept.beta, type = "diag")[[4]]

#===================================#
## Save model outputs ----
### Diagnostic plots ----
#### LDI/RLI
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-LDI-RLI-m.png"),
       cowplot::plot_grid(p.diag.lang.rli.m1, p.diag.lang.rli.m2, p.diag.lang.rli.m3, p.diag.lang.rli.m4, 
                   ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")            

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-LDI-RLI-m-int.png"),
       cowplot::plot_grid(p.diag.lang.rli.m2.1, p.diag.lang.rli.m2.2, p.diag.lang.rli.m2.3, p.diag.lang.rli.m2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")            

#### LDI ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-LDI.png"),
       cowplot::plot_grid(p.diag.lang1, p.diag.lang2, p.diag.lang3, p.diag.lang4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-LDI-int.png"),
       cowplot::plot_grid(p.diag.lang2.1, p.diag.lang2.2, p.diag.lang2.3, p.diag.lang2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

#### Animals ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-animal.png"),
       cowplot::plot_grid(p.diag.animal1, p.diag.animal2, p.diag.animal3, p.diag.animal4, 
                   ncol = 2),
width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-animal-int.png"),
       cowplot::plot_grid(p.diag.animal2.1, p.diag.animal2.2, p.diag.animal2.3, p.diag.animal2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

#### Amphibians ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-amph.png"),
       cowplot::plot_grid(p.diag.amph1, p.diag.amph2, p.diag.amph3, p.diag.amph4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-amph-int.png"),
       cowplot::plot_grid(p.diag.amph2.1, p.diag.amph2.2, p.diag.amph2.3, p.diag.amph2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

#### Birds ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-bird.png"),
       cowplot::plot_grid(p.diag.bird1, p.diag.bird2, p.diag.bird3, p.diag.bird4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-bird-int.png"),
       cowplot::plot_grid(p.diag.bird2.1, p.diag.bird2.2, p.diag.bird2.3, p.diag.bird2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

#### Mammals ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-mamm.png"),
       cowplot::plot_grid(p.diag.mamm1, p.diag.mamm2, p.diag.mamm3, p.diag.mamm4, 
                   ncol = 2),
width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-mamm-int.png"),
       cowplot::plot_grid(p.diag.mamm2.1, p.diag.mamm2.2, p.diag.mamm2.3, p.diag.mamm2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")

#### Reptiles ----
ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-rept.png"),
       cowplot::plot_grid(p.diag.rept1, p.diag.rept2, p.diag.rept3, p.diag.rept4, 
                   ncol = 2),
width = 15, height = 10, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-DIAGNOSTICS-RLI-rept-int.png"),
       cowplot::plot_grid(p.diag.rept2.1, p.diag.rept2.2, p.diag.rept2.3, p.diag.rept2.4, 
                          ncol = 2),
       width = 15, height = 10, unit = "in", bg = "white")





#### Coefficient plot for all models
## Beta regression
## Occupation time only model ----
all.models <- list()
all.models[[1]] <- m.resid.LDI.RLI.m.beta
all.models[[2]] <- m.resid.lang.beta
all.models[[3]] <- m.resid.animal.beta


p.coeff.all <- plot_models(all.models,
                           title = "A) Residual analysis (without interaction)",
                           show.values = F,
                           auto.label = T,
                           vline.color = "darkgrey",
                           p.shape = T,
                           spacing = 0.6,
                           #colors = c("#0243b3", "#199ebd", "#b21807", "#7d12e4", "#087224", "#675128", "#ee5533"),
                           colors = c("#199ebd", "#7d12e4", "#087224"),
                           m.labels = c("LDI/RLI", "LDI", "RLI"),
                           axis.lim = c(-0.5, 0.5),
                           legend.title = "",
                           dot.size	= 6,
                           line.size = 1) + 
  scale_x_discrete(labels=rev(c("Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time"))) +
  theme_bw() +
  theme(text = element_text(size = 22),
        legend.position = "") + 
  ylim(-0.7, 0.7)

## Occupation time, number colonizers interaction model ----
all.models.int <- list()
all.models.int[[1]] <- m.resid.LDI.RLI.m2.beta
all.models.int[[2]] <- m.resid.lang2.beta
all.models.int[[3]] <- m.resid.animal2.beta

p.coeff.all.int <- plot_models(all.models.int,
                               title = "B) Residual analysis (with interaction)",
                               show.values = F,
                               vline.color = "black",
                               p.shape = T,
                               spacing = 0.6,
                               colors = c("#199ebd", "#7d12e4", "#087224"),
                               m.labels = c("LDI/RLI", "LDI", "RLI"),
                               axis.lim = c(-0.5, 0.5),
                               dot.size	= 6,
                               line.size = 1) + 
  scale_x_discrete(labels=rev(c("Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time x\nnumber of colonizers"))) +
  theme_bw() +
  theme(text = element_text(size = 22),
        legend.position = "")+ 
  ylim(-0.7, 0.7)


legend_b <- get_legend(p.coeff.all + 
                         guides(color = guide_legend(nrow = 1)) +
                         theme(text = element_text(size = 22),
                               legend.position = "bottom",
                               legend.box="vertical", legend.margin=margin())
)

x11()
cowplot::plot_grid(cowplot::plot_grid(p.coeff.all, p.coeff.all.int, ncol = 2),
                   legend_b, ncol = 1, rel_heights = c(1,0.1))


## Save model output ----
ggsave(filename = paste0(path.fig, "/", datum, "-Figure-Resid-Model-coefficients.png"),
      p.coeff.all,
      width = 10, height = 9, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-Figure-Resid-Model-coefficients-int.png"),
       p.coeff.all.int,
       width = 10, height = 9, unit = "in", bg = "white")

ggsave(filename = paste0(path.fig, "/", datum, "-Figure-Resid-Model-coefficients-combined.png"),
       cowplot::plot_grid(cowplot::plot_grid(p.coeff.all, p.coeff.all.int, ncol = 2),
                   legend_b, ncol = 1, rel_heights = c(1,0.1)),
width = 18, height = 10, unit = "in", bg = "white")

#### Model summary table for all models
tab_model(m.resid.LDI.RLI.m.beta, m.resid.lang.beta, m.resid.animal.beta,
          pred.labels = c("Intercept", "Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time"),
          linebreak = T,
          dv.labels = c("LDI/RLI indicator", "LDI", "RLI"),
          title = "A)", 
          file =paste0(path.fig, "/", datum, "-Resid-Model-Ind-LDI-RLI.html"))


tab_model(m.resid.LDI.RLI.m2.beta, m.resid.lang2.beta, m.resid.animal2.beta,
          pred.labels = c("Intercept", "Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time x number of colonizers"),
          linebreak = T,
          dv.labels = c("LDI/RLI indicator", "LDI", "RLI"), 
          title = "b)", 
          file =paste0(path.fig, "/", datum, "-Resid-Model-Ind-LDI-RLI-int.html"))


tab_model(m.resid.amph.beta, m.resid.bird.beta, m.resid.mamm.beta, m.resid.rept.beta,
          pred.labels = c("Intercept", "Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time"),
          linebreak = T,
          dv.labels = c("RLI (amphibians)", "RLI (birds)", "RLI (mammals)", "RLI (reptiles)"), 
          title = "A)", 
          file =paste0(path.fig, "/", datum, "-Resid-Model-RLI-tax.html"))

tab_model(m.resid.amph2.beta, m.resid.bird2.beta, m.resid.mamm2.beta, m.resid.rept2.beta,
          pred.labels = c("Intercept", "Roughness", "Aridity Index", "Remoteness", "Prop. urban area", "Prop. intensive agriculture", "Prop. mixed agriculture", "GDPpc", "Years of schooling", "Occupation time x number of colonizers"),
          linebreak = T,
          dv.labels = c("RLI (amphibians)", "RLI (birds)", "RLI (mammals)", "RLI (reptiles)"), 
          title = "B)", 
          file =paste0(path.fig, "/", datum, "-Resid-Model-RLI-tax-int.html"))



tab_model(m.resid.LDI.RLI.m.beta, m.resid.lang.beta, m.resid.animal.beta, m.resid.amph.beta, m.resid.bird.beta, m.resid.mamm.beta, m.resid.rept.beta, 
          file =paste0(path.fig, "/", datum, "-Resid-Model-Summary.html"))


tab_model(m.resid.LDI.RLI.m2.beta, m.resid.lang2.beta, m.resid.animal2.beta, m.resid.amph2.beta, m.resid.bird2.beta, m.resid.mamm2.beta, m.resid.rept2.beta,
          file =paste0(path.fig, "/", datum, "-Resid-Model-Summary-int.html"))


#=======================================================================================#
#=======================================================================================#
#=======================================================================================#
#=======================================================================================#
# Hotspot visualization ----
## Calculate hotspot quantiles ----
### LDI/RLI (mean) hotspots ----
quant.LDI.RLI.m <- quantile(df.biocultural$Ind.m.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.biocultural$LDI.RLI.m.resid.cat <- NA
df.biocultural[which(df.biocultural$Ind.m.resid.beta <= quant.LDI.RLI.m[1]),]$LDI.RLI.m.resid.cat <- 1
df.biocultural[which(df.biocultural$Ind.m.resid.beta > quant.LDI.RLI.m[1] & df.biocultural$Ind.m.resid.beta <= quant.LDI.RLI.m[2]),]$LDI.RLI.m.resid.cat <- 2
df.biocultural[which(df.biocultural$Ind.m.resid.beta > quant.LDI.RLI.m[2] & df.biocultural$Ind.m.resid.beta <= quant.LDI.RLI.m[3]),]$LDI.RLI.m.resid.cat <- 3
df.biocultural[which(df.biocultural$Ind.m.resid.beta > quant.LDI.RLI.m[3] & df.biocultural$Ind.m.resid.beta <= quant.LDI.RLI.m[4]),]$LDI.RLI.m.resid.cat <- 4
df.biocultural[which(df.biocultural$Ind.m.resid.beta > quant.LDI.RLI.m[4] & df.biocultural$Ind.m.resid.beta <= quant.LDI.RLI.m[5]),]$LDI.RLI.m.resid.cat <- 5
df.biocultural[which(df.biocultural$Ind.m.resid.beta > quant.LDI.RLI.m[5]),]$LDI.RLI.m.resid.cat <- 6

### LDI hotspots ----
quant.LDI <- quantile(df.linguistic$LDI.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.linguistic$LDI.resid.cat <- NA
df.linguistic[which(df.linguistic$LDI.resid.beta <= quant.LDI[1]),]$LDI.resid.cat <- 1
df.linguistic[which(df.linguistic$LDI.resid.beta > quant.LDI[1] & df.linguistic$LDI.resid.beta <= quant.LDI[2]),]$LDI.resid.cat <- 2
df.linguistic[which(df.linguistic$LDI.resid.beta > quant.LDI[2] & df.linguistic$LDI.resid.beta <= quant.LDI[3]),]$LDI.resid.cat <- 3
df.linguistic[which(df.linguistic$LDI.resid.beta > quant.LDI[3] & df.linguistic$LDI.resid.beta <= quant.LDI[4]),]$LDI.resid.cat <- 4
df.linguistic[which(df.linguistic$LDI.resid.beta > quant.LDI[4] & df.linguistic$LDI.resid.beta <= quant.LDI[5]),]$LDI.resid.cat <- 5
df.linguistic[which(df.linguistic$LDI.resid.beta > quant.LDI[5]),]$LDI.resid.cat <- 6

### Animal hotspots ----
quant.animal <- quantile(df.biodiversity$RLI.animal.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.biodiversity$RLI.animal.resid.cat <- NA
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta <= quant.animal[1]),]$RLI.animal.resid.cat <- 1
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta > quant.animal[1] & df.biodiversity$RLI.animal.resid.beta <= quant.animal[2]),]$RLI.animal.resid.cat <- 2
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta > quant.animal[2] & df.biodiversity$RLI.animal.resid.beta <= quant.animal[3]),]$RLI.animal.resid.cat <- 3
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta > quant.animal[3] & df.biodiversity$RLI.animal.resid.beta <= quant.animal[4]),]$RLI.animal.resid.cat <- 4
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta > quant.animal[4] & df.biodiversity$RLI.animal.resid.beta <= quant.animal[5]),]$RLI.animal.resid.cat <- 5
df.biodiversity[which(df.biodiversity$RLI.animal.resid.beta > quant.animal[5]),]$RLI.animal.resid.cat <- 6


df.linguistic %>% filter(LDI.resid.cat == 1) %>% select(Alpha.3.code, LDI.new, LDI.resid.beta, LDI.resid.cat) %>% arrange(LDI.resid.beta)
df.biodiversity %>% filter(RLI.animal.resid.cat == 1) %>% select(Alpha.3.code, RLI.animal, RLI.animal.resid.beta, RLI.animal.resid.cat) %>% arrange(RLI.animal.resid.beta)
df.biocultural %>% filter(LDI.RLI.m.resid.cat == 1) %>% select(Alpha.3.code, Ind.mean, Ind.m.resid.beta, LDI.RLI.m.resid.cat) %>% arrange(Ind.m.resid.beta)

df.linguistic %>% filter(LDI.resid.cat == 6) %>% select(English.short.name, Alpha.3.code, LDI.new, LDI.resid.beta, LDI.resid.cat) %>% arrange(LDI.resid.beta)
df.biodiversity %>% filter(RLI.animal.resid.cat == 6) %>% select(English.short.name, Alpha.3.code, RLI.animal, RLI.animal.resid.beta, RLI.animal.resid.cat) %>% arrange(desc(RLI.animal.resid.beta))
df.biocultural %>% filter(LDI.RLI.m.resid.cat == 6) %>% select(English.short.name, Alpha.3.code, Ind.mean, Ind.m.resid.beta, LDI.RLI.m.resid.cat) %>% arrange(desc(Ind.m.resid.beta))


### Amphibian hotspots ----
quant.amph <- quantile(df.rli.amph.s$RLI.amph.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.rli.amph.s$RLI.amph.resid.cat <- NA
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta <= quant.amph[1]),]$RLI.amph.resid.cat <- 1
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta > quant.amph[1] & df.rli.amph.s$RLI.amph.resid.beta <= quant.amph[2]),]$RLI.amph.resid.cat <- 2
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta > quant.amph[2] & df.rli.amph.s$RLI.amph.resid.beta <= quant.amph[3]),]$RLI.amph.resid.cat <- 3
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta > quant.amph[3] & df.rli.amph.s$RLI.amph.resid.beta <= quant.amph[4]),]$RLI.amph.resid.cat <- 4
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta > quant.amph[4] & df.rli.amph.s$RLI.amph.resid.beta <= quant.amph[5]),]$RLI.amph.resid.cat <- 5
df.rli.amph.s[which(df.rli.amph.s$RLI.amph.resid.beta > quant.amph[5]),]$RLI.amph.resid.cat <- 6

### Bird hotspots ----
quant.bird <- quantile(df.rli.bird.s$RLI.bird.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.rli.bird.s$RLI.bird.resid.cat <- NA
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta <= quant.bird[1]),]$RLI.bird.resid.cat <- 1
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta > quant.bird[1] & df.rli.bird.s$RLI.bird.resid.beta <= quant.bird[2]),]$RLI.bird.resid.cat <- 2
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta > quant.bird[2] & df.rli.bird.s$RLI.bird.resid.beta <= quant.bird[3]),]$RLI.bird.resid.cat <- 3
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta > quant.bird[3] & df.rli.bird.s$RLI.bird.resid.beta <= quant.bird[4]),]$RLI.bird.resid.cat <- 4
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta > quant.bird[4] & df.rli.bird.s$RLI.bird.resid.beta <= quant.bird[5]),]$RLI.bird.resid.cat <- 5
df.rli.bird.s[which(df.rli.bird.s$RLI.bird.resid.beta > quant.bird[5]),]$RLI.bird.resid.cat <- 6

### Mammal hotspots ----
quant.mamm <- quantile(df.rli.mamm.s$RLI.mamm.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.rli.mamm.s$RLI.mamm.resid.cat <- NA
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta <= quant.mamm[1]),]$RLI.mamm.resid.cat <- 1
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta > quant.mamm[1] & df.rli.mamm.s$RLI.mamm.resid.beta <= quant.mamm[2]),]$RLI.mamm.resid.cat <- 2
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta > quant.mamm[2] & df.rli.mamm.s$RLI.mamm.resid.beta <= quant.mamm[3]),]$RLI.mamm.resid.cat <- 3
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta > quant.mamm[3] & df.rli.mamm.s$RLI.mamm.resid.beta <= quant.mamm[4]),]$RLI.mamm.resid.cat <- 4
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta > quant.mamm[4] & df.rli.mamm.s$RLI.mamm.resid.beta <= quant.mamm[5]),]$RLI.mamm.resid.cat <- 5
df.rli.mamm.s[which(df.rli.mamm.s$RLI.mamm.resid.beta > quant.mamm[5]),]$RLI.mamm.resid.cat <- 6

### Reptile hotspots ----
quant.rept <- quantile(df.rli.rept.s$RLI.rept.resid.beta, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

df.rli.rept.s$RLI.rept.resid.cat <- NA
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta <= quant.rept[1]),]$RLI.rept.resid.cat <- 1
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta > quant.rept[1] & df.rli.rept.s$RLI.rept.resid.beta <= quant.rept[2]),]$RLI.rept.resid.cat <- 2
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta > quant.rept[2] & df.rli.rept.s$RLI.rept.resid.beta <= quant.rept[3]),]$RLI.rept.resid.cat <- 3
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta > quant.rept[3] & df.rli.rept.s$RLI.rept.resid.beta <= quant.rept[4]),]$RLI.rept.resid.cat <- 4
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta > quant.rept[4] & df.rli.rept.s$RLI.rept.resid.beta <= quant.rept[5]),]$RLI.rept.resid.cat <- 5
df.rli.rept.s[which(df.rli.rept.s$RLI.rept.resid.beta > quant.rept[5]),]$RLI.rept.resid.cat <- 6




#==============================================================================#
# Summarize descriptive info ----
## Hotspots & coldspots LDI/RLI indicator ----
df.biocultural.s %>% 
  filter(LDI.RLI.m.resid.cat %in% c(1,6)) %>%
  select(region,English.short.name, Alpha.3.code, Ind.mean, Ind.m.resid.beta, LDI.RLI.m.resid.cat) %>%
  arrange(Ind.m.resid.beta)

## Hotspots & coldspots LDI  
df.linguistic.s %>% 
  filter(LDI.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, LDI, LDI.resid.beta, LDI.resid.cat) %>%
  arrange(LDI.resid.beta)

## Hotspots & coldspots RLI  
df.biodiversity.s %>% 
  filter(RLI.animal.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, RLI.animal, RLI.animal.resid.beta, RLI.animal.resid.cat) %>%
  arrange(RLI.animal.resid.beta)


## Hotspots & coldspots RLI amphibians  
df.rli.amph.s %>% 
  filter(RLI.amph.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, RLI.amph, RLI.amph.resid.beta, RLI.amph.resid.cat) %>%
  arrange(RLI.amph.resid.beta)
df.rli.bird.s %>% 
  filter(RLI.bird.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, RLI.birds, RLI.bird.resid.beta, RLI.bird.resid.cat) %>%
  arrange(RLI.bird.resid.beta)
df.rli.mamm.s %>% 
  filter(RLI.mamm.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, RLI.mamm, RLI.mamm.resid.beta, RLI.mamm.resid.cat) %>%
  arrange(RLI.mamm.resid.beta)
df.rli.rept.s %>% 
  filter(RLI.rept.resid.cat %in% c(1,6)) %>%
  select(region, English.short.name, Alpha.3.code, RLI.rept, RLI.rept.resid.beta, RLI.rept.resid.cat) %>%
  arrange(RLI.rept.resid.beta)


#==============================================================================#
# Calculate differences between corrected LDI and RLI ----
df.LDI.RLI <- df.biocultural.s %>%
  left_join(df.linguistic.s %>% select(Alpha.3.code, LDI.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.biodiversity.s %>% select(Alpha.3.code, RLI.animal.resid.beta), by = "Alpha.3.code")


quantile(df.LDI.RLI$LDI)
quantile(df.LDI.RLI$RLI.animal)

x11()
boxplot(df.LDI.RLI$LDI, df.LDI.RLI$RLI.animal)
boxplot(df.LDI.RLI$LDI.resid, df.LDI.RLI$RLI.animal.resid)







#####################################
# Get references for R and all packages used ----
library("grateful")
cite_packages(pkgs = "Session", 
              out.format = "docx", #citation.style = "nature",
              out.dir = path.fig)




save.image(paste0(path.data, "/", datum, "-Models-Workspace.RData")) # creating ".RData"

