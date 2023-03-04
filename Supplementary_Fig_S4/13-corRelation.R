# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##########################################################################

library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)

##########################################################################

rawdf = read.table('12-complexity+stability.txt', header = T, sep = '\t', 
                   stringsAsFactors = F, row.names = 1)
rawdf.c = t(rawdf[, str_detect(colnames(rawdf), 'C')])
rawdf.m = t(rawdf[, str_detect(colnames(rawdf), 'M')])

cc = cor(rawdf.c)
cm = cor(rawdf.m)

cc = as.data.frame(cc[1:8, 9:11])
cm = as.data.frame(cm[1:8, 9:11])

cc$complex = rownames(cc)
cm$complex = rownames(cm)

cc.r = melt(cc, variable.name = "Stability", measure.name = colnames(cc)[1:3],  
            value.name ="cor")
cc.r$cor = round(cc.r$cor, 2)
cm.r = melt(cm, variable.name = "Stability", measure.name = colnames(cm)[1:3],  
            value.name ="cor")
cm.r$cor = round(cm.r$cor, 2)

cc.r$complex = factor(cc.r$complex, levels = rev(rownames(cc)))
cm.r$complex = factor(cm.r$complex, levels = rev(rownames(cm)))

cc.r$text = cc.r$cor
cc.r$text[abs(cc.r$text) < 0.7] = NA
pc = ggplot(cc.r, aes(Stability, complex, fill = cor))+
  geom_raster(color = "white")+
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",  
                       midpoint = 0, limit = c(-1,1), space = "Lab") +
  geom_text(aes(Stability, complex, label = text), color = "black", size = 2) +
  theme( panel.background = element_blank(),
         panel.grid = element_blank(), 
         axis.title = element_blank(), 
         axis.ticks = element_blank(), 
         axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1, color = 'black'),
         axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, color = 'black'),
         legend.key.width = unit(.07, 'in'),
         legend.key.height = unit(.2, 'in'),
         legend.position = "right",
         legend.background = element_blank(),
         legend.key = element_blank())

cm.r$text = cm.r$cor
cm.r$text[abs(cm.r$text) < 0.7] = NA
pm = ggplot(cm.r, aes(Stability, complex, fill = cor))+
  geom_raster(color = "white")+
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",  
                       midpoint = 0, limit = c(-1,1), space = "Lab") +
  geom_text(aes(Stability, complex, label = text), color = "black", size = 2) +
  theme( panel.background = element_blank(),
         panel.grid = element_blank(), 
         axis.title = element_blank(), 
         axis.ticks = element_blank(), 
         axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1, color = 'black'),
         axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, color = 'black'),
         legend.key.width = unit(.07, 'in'),
         legend.key.height = unit(.2, 'in'),
         legend.position = "right",
         legend.background = element_blank(),
         legend.key = element_blank())


p = grid.arrange(pc, pm, nrow = 1)

ggsave('14-cor_complex_stability.pdf',  p, units = 'in', width = 5, height = 2.5)

