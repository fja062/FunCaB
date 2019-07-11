#Package to make plots
library(RColorBrewer)
library(cowplot)

#source("/Users/fja062/Documents/seedclimComm/seedclimComm/inst/graminoidRemovals/multiplot_function.R")


#### palettes and labelling ####
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#1C9099", "#A6BDDB", "#ECE2F0", "orange3", "white")

plot(1:length(cbPalette), col = cbPalette, pch = 16, cex = 5) #check colour and order

precip.lab <-   scale_x_discrete("Precipitation [mm y-1]",
                                 labels = c("0.6"="600", "1.2"="1200","2"="2000","2.7"="2700"))
temp.lab <- scale_x_discrete("Temperature [C]",
                             labels = c("0.6"="600", "1.2"="1200","2"="2000","2.7"="2700"))

axis.dim <- theme(axis.text=element_text(size=11),
                  axis.title=element_text(size=14),
                  axis.ticks = element_blank(),
                  legend.text = element_text(size=11),
                  legend.title = element_text(size=13),
                  strip.text.x = element_text(size = 13),
                  strip.text.y = element_text(size = 13))

axis.dimLarge <- theme(axis.text=element_text(size=13),
                  axis.title=element_text(size=15),
                  axis.ticks = element_blank(),
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=15),
                  strip.text.x = element_text(size = 15),
                  strip.text.y = element_text(size = 15))


legend.title.prec <- "Treatment and \n precipitation"
legend.title.temp <- "Treatment \n and temperature"
legend.title.treat <- "Treatment"
legend.title.weat <- " \n "
legend.title.climate <- ""
legend.TOD <- "Time of day"

treatOrder <- c(c("C", "FB", "GF", "GB", "FGB"))
tempLabs <- c(`6.5` = "6.5ºC", `8.5` = "8.5ºC", `10.5` = "10.5ºC")
precipLabs <- c(`0.6` = "0.6mm/1000", `1.2` = "1.2mm/1000", `2` = "2mm/1000", `2.7` = "2.7mm/1000")
