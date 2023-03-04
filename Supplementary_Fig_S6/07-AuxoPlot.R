# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##############################################################################

library(tidyverse)
library(reshape2)
library(grid)

otumodu = read.table('07-OtuModuleFilt.txt', header = T, sep = '\t', row.names = 1)

c2m = read.table('08-AuxotrophiesList.txt', header = F, sep = '\t')
rownames(c2m) = c2m$V4
setcore = intersect(c2m$V4, colnames(otumodu))
c2m = c2m[setcore, ]

newotumodu = data.frame(OTUID = rownames(otumodu))
rownames(newotumodu) = newotumodu$OTUID

for (c in unique(c2m$V3)) {
  newotumodu$new = c()
  ms = c2m$V4[c2m$V3 == c]
  for (rotu in 1:nrow(otumodu)) {
    newotumodu$new[rotu] = max(as.numeric(otumodu[rotu, ms]))
  }
  colnames(newotumodu)[ncol(newotumodu)] = c
}

write.table(newotumodu, '10-maxC2M.txt', append = F,quote = F,sep = '\t', 
            row.names = F, col.names = T)

##########################################################################

taxdf = read.table('09-GTDBtax.txt', header = T, sep = '\t')

taxdf.f = taxdf[, c('OTUID', 'Family')] 

newotumodu2 = merge(newotumodu, taxdf.f, by = 'OTUID', all.x = T) 

##########################################################################

filt.fam = as.data.frame(table(newotumodu2$Family))
fams.filt = as.character(filt.fam$Var1[filt.fam$Freq >= 20 & filt.fam$Freq <= 300 ])

sort(unique(newotumodu2$Family))

fams = fams.filt[!(fams.filt %in% c('Desulfobacteraceae', 'Desulfobulbaceae', 
                                    'Desulfovibrionaceae', 'Marinifilaceae', 
                                    'Marinilabiliaceae', 'Prolixibacteraceae'))]

# fams = c('Desulfobacteraceae', 'Desulfobulbaceae', 'Desulfovibrionaceae', 
#          'Marinifilaceae', 'Marinilabiliaceae', 'Prolixibacteraceae')

newotumodu3 = newotumodu2[newotumodu2$Family %in% fams, ]

write.table(data.frame(table(newotumodu3$Family)), '12-familySummary.txt', 
            append = F,quote = F,sep = '\t', 
            row.names = F, col.names = T)

##########################################################################

omdf = melt(newotumodu3, measure.vars = colnames(newotumodu3)[2:(ncol(newotumodu3)-1)])
omdf$variable = as.vector(omdf$variable)

auxolist = read.table('11-AuxoFiltOrder.txt', header = F, sep = '\t', row.names = 2)
for (auxo in rownames(auxolist)) {
  omdf$variable[omdf$variable == auxo] = as.character(auxolist[auxo,])
}

omdf$variable = factor(omdf$variable, levels = auxolist$V1)

##########################################################################

pauxo1 = ggplot(omdf, aes(x = OTUID, y = variable, fill = value)) +
  geom_tile(colour = NA) + 
  facet_grid(.~Family, scales = 'free_x', space = 'free_x') + 
  theme_bw()+
  theme(plot.margin=unit(c(0, 0, 0, 3),'mm'),
        panel.background = element_rect(fill="white", color="black", size = 0),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = 'black', hjust = 1, vjust = 0.5), 
        strip.background = element_blank(),
        strip.text.x  = element_text(size = 8, colour = 'black', angle = 90, 
                                     hjust = 0, vjust = 0.5),
        legend.key.height = unit(.05, 'in'),
        legend.key.width = unit(0.1, 'in'),
        legend.title = element_blank(),
        legend.text =  element_text(size = 6, colour = 'black'), 
        legend.position = 'bottom' )


ggsave('13-plotAuxo.pdf', pauxo1, 
       width = 6, height = 4.3, units = 'in')
