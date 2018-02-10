#########################################
#  	Helper functions	    	#
#########################################

##  Written: Sara Ballouz

## Necessary libraries
library(shiny) 
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggvis)
library(gplots)
library(pheatmap)
library(shinyBS)
library(plotly)
library(markdown)
library(NMF)
library(heatmaply)
library(readr)
library(MASS)
library(Matrix)
library(scales)
library(zoo)
library(plyr)
library(RColorBrewer)
library(ape)
#library(WGCNA)
library(dynamicTreeCut)
library(viridisLite)
library(beanplot)
library(beeswarm)
library(corrgram)

cols  = colorpanel(16, "red", "blue")
cols2 = brewer.pal(8, "Spectral")
cols3 = rainbow(30)
cols4 = colorpanel(65, "grey", "blue", "darkblue")
cols5 = colorpanel(300, "grey", "red", "darkred")
cols6 = colorpanel(30, "grey", "red", "darkmagenta")
cols7 = c("seagreen", "black", "darkmagenta")
cols8 = viridis(10)

