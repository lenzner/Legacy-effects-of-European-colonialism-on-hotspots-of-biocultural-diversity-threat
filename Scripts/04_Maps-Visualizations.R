#=======================================#
#                                       #
# Visuliazation                         #
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
packages <- c("stringr","ncdf4","tictoc","tidyverse", "sf","data.table","GGally","lme4","glmmTMB","cowplot","sjPlot","tinytable", "corrplot", "ggcorrplot", "ggpubr", "DHARMa", "gamlss", "performance", "patchwork", "grid", "gridExtra", "ggpubr") 

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Load packages
invisible(lapply(packages, library, character.only = TRUE))

#===========================================#
# Directory and paths ----
drive <- "F"
path.datasets <- paste0(drive, ':/Datasets') # path for directory with predictor datasets
path.proj <- paste0(drive, ':/Publication_Projects/2022_Linguistics_Biodiv_Language_Threat') # path for directory with code
path.data <- paste0(path.proj, '/Data') # path for directory with analysis datasets
path.fig <- paste0(path.proj, '/Figures/2024-09-Final-figures') # path for directory with analysis datasets
path.code <- paste0(path.proj, '/R-Scripts') # path for directory with analysis datasets

# Load workspace ----
load(paste0(path.data, "/2024-09-03-Hotspot-Models.RData")) # Load hotspot models
load(paste0(path.data, "/2024-09-03-Dataframes.RData")) # Load dataframes for each threat indicator
load(paste0(path.data, "/2024-04-22_Full-Analysis-Dataset.Rdata"))
load(paste0(path.data, "/2024-07-09-RedList-cat-ALL-by-region.RData")) # RedList categories by group and region


names.hotspot.models <- ls()[grep("hotspot",ls())]
names.dataframes <- ls()[grep("df",ls())]

# Functions ----
source(paste0(path.code,"/F02-Hotspot-Visualization.R"))
source(paste0(path.code,"/F03-1b1-distance.R"))

# Shapefile ----
## Final region selection (countries GADM level 0) ----
### #gadm.data.path <- paste0(path.datasets, "/Regions/GADM/gadm_410-levels.gpkg")
### gadm0 <- st_read(gadm.data.path, layer = "ADM_0") # national

gadm.data.path <- paste0(path.datasets, "/Regions/GADM/GADM0-simpl-0.1.shp")
gadm0 <- st_read(gadm.data.path) # national simplified
shape <- gadm0 %>%
  st_wrap_dateline() %>% # avoid banding
  st_transform(crs = "+proj=moll")

sf::sf_use_s2(FALSE) # disable s2 proccessing to run st_union() function
shape.proj <- shape



# Map visualizations ----
x11()
## Plot with hot- and coldspot extremes as inset
plot.hotspots(indicator = "biocultural", map = "hotspots") +
  inset_element(plot.hotspots(indicator = "biocultural", map = "extremes"), left = 0, bottom = -0.1, right = 0.35, top = 0.4)



p.biocultural.hot <- plot.hotspots(indicator = "biocultural", map = "hotspots", composite = T)
p.linguistic.hot <- plot.hotspots(indicator = "linguistic", map = "hotspots", composite = T)
p.biodiversity.hot <- plot.hotspots(indicator = "biodiversity", map = "hotspots", composite = T)

p.amph.hot <- plot.hotspots(indicator = "amph", map = "hotspots", composite = T)
p.bird.hot <- plot.hotspots(indicator = "bird", map = "hotspots", composite = T)
p.mamm.hot <- plot.hotspots(indicator = "mamm", map = "hotspots", composite = T)
p.rept.hot <- plot.hotspots(indicator = "rept", map = "hotspots", composite = T)

ggsave(filename = paste0(path.fig, "/2024-02-02-Results/", datum, "-Hotspots-LDI-RLI-combined1.pdf"),
       cowplot::plot_grid(p.linguistic.hot, p.biodiversity.hot, p.biocultural.hot, nrow = 3,
                          align = "v"),
       width = 210, height = 297, unit = "mm", bg = "white")




# Plot with only one legend
ggsave(filename = paste0(path.fig, "/2024-02-02-Results/", datum, "-Hotspots-LDI-RLI-combined.pdf"),
       grid_arrange_shared_legend(p.linguistic.hot, p.biodiversity.hot, p.biocultural.hot, nrow = 3, ncol = 1),
       width = 210, height = 297, unit = "mm", bg = "white")

  ggsave(filename = paste0(path.fig, "/2024-02-02-Results/", datum, "-Hotspots-animal-groups-combined.png"),
       grid_arrange_shared_legend(p.amph.hot, p.bird.hot, p.mamm.hot, p.rept.hot, nrow = 2, ncol = 2),
       width = 210, height = 210, unit = "mm", bg = "white")






#==============================================================================#
#==============================================================================#
## Spatial correlation pattern between linguistic and biological ----
## Drop correlated predictors ----
data.global <- data.full %>%
  select(!c(dem.range, aridity.sd, crop.area.sum, pop.2015, struct_div, pop.dens.2015,crop.dens.2019, hdi.pc.2015.mean, traveltime.mean1, traveltime.mean9,
            Alpha.2.code))

data.corr.pattern <- data.global %>%
  select(English.short.name:Numeric) %>%
  left_join(df.biocultural %>% select(Alpha.3.code, Ind.m.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.linguistic %>% select(Alpha.3.code, LDI.resid.beta), by = "Alpha.3.code") %>%
  left_join(df.biodiversity %>% select(Alpha.3.code, RLI.animal.resid.beta), by = "Alpha.3.code")


### Calculate distance to 1:1 line ----
#### LDI-RLI ----
diff.pattern.RLI.LDI <- data.corr.pattern %>% 
  mutate(dist.RLI.LDI = apply(data.corr.pattern[,c("RLI.animal.resid.beta", "LDI.resid.beta")],
                              1, function(x) dist_point_line(a = x, slope = 1, intercept = 0)))

#### Scatterplot correlation linguistic - biological diversity ----
p.diff.points <- ggplot(diff.pattern.RLI.LDI, aes(x = RLI.animal.resid.beta, y = LDI.resid.beta, colour = dist.RLI.LDI), show.legend = F) +
  geom_point(show.legend = F) +
  geom_abline (slope=1, linetype = "dashed", color="darkgrey", linewidth = 0.5) +
  scale_color_gradient2(low = "#9831CC", mid = "grey", high = "#00A693", na.value = "white") +
  ylim(-3.5,3.5) + xlim(-4.5,2) +
  ggtitle("B)") +
  ylab("RLI languages\n(residuals)") + xlab("RLI animals\n(residuals)") +
  ggrepel::geom_text_repel(aes(label = Alpha.3.code),
                           min.segment.length = 0,
                           #box.padding   = 0, 
                           #point.padding = 0,
                           #segment.color = ,
                           size = 5) +
  theme_minimal() +
  font("title", size = 14, face = "bold") +
  font("ylab", size = 12) +
  font("xlab", size = 12) +
  theme(text = element_text(size=9),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position="none")



#### Map correlation linguistic - biological diversity ----
### LDI/RLI - difference ----
shape.dist <- shape.proj %>%
  left_join(diff.pattern.RLI.LDI, by = c("GID_0" = "Alpha.3.code")) %>%
  mutate(
    area = st_area(shape.proj),
    area = as.numeric(area),
    area = round(area*0.000001,0)
  )

shape.dist.isl <- shape.dist %>%
  filter(area < 35000) %>%
  st_centroid()


p.diff.hot <- ggplot() +
  geom_sf(data = shape.dist, aes(geometry = geometry, group = dist.RLI.LDI, fill = dist.RLI.LDI), show.legend = T, color = "lightgrey", size = 0.05) +
  scale_fill_gradient2(low = "#9831CC", mid = "grey", high = "#00A693", na.value = "white", name = "1:1 distance") +
  geom_sf(data = shape.dist.isl, aes(geometry = geometry, group = dist.RLI.LDI, color = dist.RLI.LDI), size = 3, shape = "O", show.legend = FALSE) +
  #scale_fill_continuous(type = "gradient", mid = "white", low = "#943126", high = "#2874A6", name = "1:1 distance") +
  theme_minimal() + 
  theme(text = element_text(size=10)) +
  labs(title = "A)") +
  font("title", size = 14, face = "bold") +
  font("legend.title", size = 12) +
  font("legend.text", size = 12)


cowplot::plot_grid(p.diff.points, p.diff.hot, ncol = 1)



ggsave(filename = paste0(path.fig, "/", datum, "-CORRELATION-linguistic-biodiversity-hotspots.pdf"),
       cowplot::plot_grid(p.diff.points, p.diff.hot, ncol = 1, align = "v"),
       width = 210, height = 198, unit = "mm", bg = "white")





















#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
### Insets ----
#### Inset Europe ----
target_crs <- '+proj=moll'

# Coordinates in WGS coordinates
disp_win_wgs84 <- st_sfc(st_point(c(-20, 30)), st_point(c(45, 73)),
                         crs = 4326)
disp_win_wgs84

# Transform coordinates to Mollweide crs
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_trans

# Get coordinates to pass as limits for inset
disp_win_coord <- st_coordinates(disp_win_trans)


p.hot.LDI.RLI.inset.EU <- ggplot() + geom_sf(data = shape.lang.RLI.hot, aes(geometry = geom, group = LDI.RLI.m.resid.cat, fill = factor(LDI.RLI.m.resid.cat)), color = "lightgrey", size = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = cols, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.lang.RLI.hot.isl, aes(geometry = geom, group = LDI.RLI.m.resid.cat, size = 3, color = as.factor(LDI.RLI.m.resid.cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size=10)) +
  labs(title = "D) Europe")


#### Inset Central America ----
# Coordinates in WGS coordinates
disp_win_wgs84 <- st_sfc(st_point(c(-120, 5)), st_point(c(-50, 35)),
                         crs = 4326)
disp_win_wgs84

# Transform coordinates to Mollweide crs
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_trans

# Get coordinates to pass as limits for inset
disp_win_coord <- st_coordinates(disp_win_trans)


p.hot.LDI.RLI.inset.CA <- ggplot() + geom_sf(data = shape.lang.RLI.hot, aes(geometry = geom, group = LDI.RLI.m.resid.cat, fill = factor(LDI.RLI.m.resid.cat)), color = "lightgrey", size = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = cols, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.lang.RLI.hot.isl, aes(geometry = geom, group = LDI.RLI.m.resid.cat, size = 3, color = as.factor(LDI.RLI.m.resid.cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size=10)) +
  labs(title = "B) Central America")


#### Inset Indo-Maly ----
# Coordinates in WGS coordinates
disp_win_wgs84 <- st_sfc(st_point(c(75, -25)), st_point(c(180, 28)),
                         crs = 4326)
disp_win_wgs84

# Transform coordinates to Mollweide crs
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_trans

# Get coordinates to pass as limits for inset
disp_win_coord <- st_coordinates(disp_win_trans)


p.hot.LDI.RLI.inset.IM <- ggplot() + geom_sf(data = shape.lang.RLI.hot, aes(geometry = geom, group = LDI.RLI.m.resid.cat, fill = factor(LDI.RLI.m.resid.cat)), color = "lightgrey", size = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = cols, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.lang.RLI.hot.isl, aes(geometry = geom, group = LDI.RLI.m.resid.cat, size = 3, color = as.factor(LDI.RLI.m.resid.cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size=10)) +
  labs(title = "E) Malay Archipelago")

#### Southern Africa ----
# Coordinates in WGS coordinates
disp_win_wgs84 <- st_sfc(st_point(c(0, -40)), st_point(c(60, 0)),
                         crs = 4326)
disp_win_wgs84

# Transform coordinates to Mollweide crs
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_trans

# Get coordinates to pass as limits for inset
disp_win_coord <- st_coordinates(disp_win_trans)


p.hot.LDI.RLI.inset.SA <- ggplot() + geom_sf(data = shape.lang.RLI.hot, aes(geometry = geom, group = LDI.RLI.m.resid.cat, fill = factor(LDI.RLI.m.resid.cat)), color = "lightgrey", size = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = cols, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.lang.RLI.hot.isl, aes(geometry = geom, group = LDI.RLI.m.resid.cat, size = 3, color = as.factor(LDI.RLI.m.resid.cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size=10)) +
  labs(title = "C) Southern Africa")

