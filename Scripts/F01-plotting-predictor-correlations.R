#========================================================================#
# Publication Biocultural Threat                                         #
# Function for plotting predictor correlations across different subsets  #
# Bernd Lenzner                                                          #
# bernd.lenzner@univie.ac.at                                             #
#========================================================================#

plot.predictor.correlation <- function(indicator = c("lang.RLI", "LDI", "RLI.animal", "rli.amph", "rli.bird", "rli.mamm", "rli.rept", "rli.no.amph", "rli.no.bird", "rli.no.mamm", "rli.no.rept")){

  if(indicator == "lang.RLI"){dat <- df.biocultural}
  if(indicator == "lang.RLI"){title = "A) Biocultural Diversity"}
  
  if(indicator == "LDI"){dat <- df.linguistic}
  if(indicator == "LDI"){title = "B) Linguistic Diversity"}
  
  if(indicator == "RLI.animal"){dat <- df.biodiversity}
  if(indicator == "RLI.animal"){title = "C) Biological Diversity"}
  
  if(indicator == "RLI.amph"){dat <- df.rli.amph}
  if(indicator == "RLI.amph"){title = "A) Amphibians"}
  
  if(indicator == "RLI.bird"){dat <- df.rli.bird}
  if(indicator == "RLI.bird"){title = "B) Birds"}
  
  if(indicator == "RLI.mamm"){dat <- df.rli.mamm}
  if(indicator == "RLI.mamm"){title = "C) Mammals"}
  
  if(indicator == "RLI.rept"){dat <- df.rli.rept}
  if(indicator == "RLI.rept"){title = "D) Reptiles"}
  
  if(indicator == "RLI.no.amph"){dat <- df.rli.no.amph}
  if(indicator == "RLI.no.amph"){title = "Red List Index (without amphibians)"}
  
  if(indicator == "RLI.no.bird"){dat <- df.rli.no.bird}
  if(indicator == "RLI.no.bird"){title = "Red List Index (wihthout birds)"}
  
  if(indicator == "RLI.no.mamm"){dat <- df.rli.no.mamm}
  if(indicator == "RLI.no.mamm"){title = "Red List Index (without mammals)"}
  
  if(indicator == "RLI.no.rept"){dat <- df.rli.no.rept}
  if(indicator == "RLI.no.rept"){title = "Red List Index (without reptiles)"}
  
  dat <- dat %>%
    rename("Sampling bias" = samp.bias.s,
           "Prop. intense agriculture" = intense.s,
           "Roughness" = dem.s,
           "Area" = area.s,
           "Years of schooling" = schooling.s,
           "GDPpc" = gdp.s,
           "Aridity Index" = aridity.s,
           "Prop. urban area" = urban.s,
           "Occupation time" = occ.time.full.s,
           "Temperature" = temp.s,
           "Traveltime (small)" = traveltime.small.s,
           "Traveltime (large)" = traveltime.large.s,
           "Prop. mixed agriculture" = mixed.s)
  
  p <- ggcorrplot(round(cor(dat %>%
                         select(Area : `Occupation time`)),1),
           hc.order = TRUE, type = "lower",
           title = title,
           outline.col = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"),
           lab = TRUE,
           lab_size = 2, # values in boxes
           tl.cex = 9 # tick labels
           ) #+
  #font("title", size = 16, face = "bold") +
  #font("legend.title", size = 14) +
  #font("legend.text", size = 14)

  return(p)

}



