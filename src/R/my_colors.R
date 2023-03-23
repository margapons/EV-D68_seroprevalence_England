
library(RColorBrewer)
library(wesanderson)
library(rcartocolor)

cutoff_cols = brewer.pal(n=9,"Set1")[c(2,1)]

years_cols = wes_palette(n=5, name="Darjeeling1",type="discrete")[c(1,2,5,3)]

prog_2006_2017 = colorRampPalette(rcartocolor::carto_pal(name="ag_Sunset"))(12)
