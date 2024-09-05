#========================================================================#
# Publication Biocultural Threat                                         #
# Function for plotting predictor correlations across different subsets  #
# Bernd Lenzner                                                          #
# bernd.lenzner@univie.ac.at                                             #
#========================================================================#


plot.hotspots <- function(indicator = c("biocultural", "linguistic", "biodiversity", "amph", "bird", "mamm", "rept"), map = c("hotspots","extremes"), composite = FALSE){

  if(indicator == "biocultural" & composite == F){title = "Biocultural diversity"}
  if(indicator == "biocultural" & composite == T){title = "C) Biocultural diversity"}
  
  if(indicator == "linguistic" & composite == F){title = "Linguistic diversity"}
  if(indicator == "linguistic" & composite == T){title = "A) Linguistic diversity"}
  
  if(indicator == "biodiversity" & composite == F){title = "Biological diversity"}
  if(indicator == "biodiversity" & composite == T){title = "B) Biological diversity"}
  
  if(indicator == "amph"){title = "A) Amphibians"}
  
  if(indicator == "bird"){title = "B) Birds"}
  
  if(indicator == "mamm"){title = "C) Mammals"}
  
  if(indicator == "rept"){title = "D) Reptiles"}
  
  
d <- names.dataframes[grep(indicator, names.dataframes)]

dat <- get(d)

quant <- quantile(unlist(dat[grep("resid", names(dat), value=T)]), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

dat$cat <- NA
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) <= quant[1]),]$cat <- 1
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) > quant[1] & unlist(dat[grep("resid", names(dat), value=T)]) <= quant[2]),]$cat <- 2
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) > quant[2] & unlist(dat[grep("resid", names(dat), value=T)]) <= quant[3]),]$cat <- 3
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) > quant[3] & unlist(dat[grep("resid", names(dat), value=T)]) <= quant[4]),]$cat <- 4
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) > quant[4] & unlist(dat[grep("resid", names(dat), value=T)]) <= quant[5]),]$cat <- 5
dat[which(unlist(dat[grep("resid", names(dat), value=T)]) > quant[5]),]$cat <- 6


shape.hot <- shape.proj %>%
  left_join(dat, by = c("GID_0" = "Alpha.3.code"))

shape.hot.isl <- shape.hot %>%
  filter(area < 35000) %>%
  st_centroid()


#=======================================================================================#
cols <- c("#943126","#E74C3C", "#F5B7B1", "#D1F2EB","#AED6F1","#2874A6" )
cols.extr <- c("#943126","lightgrey", "lightgrey", "lightgrey","lightgrey","#2874A6" )
legdtxt <- c("lower 2.5%", "lower 10%", "lower 50%", "upper 50%", "upper 10%", "upper 2.5%")


## Hotspot maps ----
p.hot.extremes <- ggplot() +
  geom_sf(data = shape.hot, aes(geometry = geometry, group = cat, fill = factor(cat)), color = "lightgrey", size = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = cols.extr, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.hot.isl, aes(geometry = geometry, group = cat, size = 3, color = as.factor(cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols.extr) +
  theme_minimal() +
  theme(legend.position = "",
        text = element_text(size=11)) +
  guides(fill = guide_legend(nrow = 1,
                             label.position = "")) +
  labs(title = "Hot & Coldspots") +
  font("title", size = 14) +
  font("legend.title", size = 12) +
  font("legend.text", size = 12)


p.hot <- ggplot() +
  geom_sf(data = shape.hot, aes(geometry = geometry, group = cat, fill = factor(cat)), color = "lightgrey", size = 0.05) +
  scale_fill_manual(values = cols, na.value = "lightgrey", 
                    name = "Quantiles", labels = legdtxt)  +
  geom_sf(data = shape.hot.isl, aes(geometry = geometry, group = cat, size = 3, color = as.factor(cat)), shape = "O", show.legend = FALSE) +
  scale_color_manual(values = cols) +
  #geom_rect(aes(xmin = -1833617, xmax = 2065900, ymin = 3643854, ymax = 8018072), color = "darkgrey", lty = 2, fill = NA)  + # Rectangle for Europe
  #geom_rect(aes(xmin = 0, xmax = 6013365, ymin = -4789399, ymax = 0), color = "darkgrey", lty = 2, fill = NA)  + # Rectangle for Southern Africa
  #geom_rect(aes(xmin = 7072957, xmax = 16702200, ymin = -3053323, ymax = 3408870), color = "darkgrey", lty = 2, fill = NA)  + # Rectangle for Indo Malay
  #geom_rect(aes(xmin = -11998476, xmax = -4427943, ymin = 617923.6, ymax = 4223222.6), color = "darkgrey", lty = 2, fill = NA)  + # Rectangle for Central America
  theme_minimal() + 
  theme(legend.position = "bottom",
        text = element_text(size=11)) +
  guides(fill = guide_legend(nrow = 1,
                             label.position = "bottom")) +
  labs(title = title) +
  font("title", size = 16, face = "bold") +
  font("legend.title", size = 12) +
  font("legend.text", size = 12)

if(map == "hotspots") return(p.hot)
if(map == "extremes") return(p.hot.extremes)

}



grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

