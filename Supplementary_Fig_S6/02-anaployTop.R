# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##############################################################################

library(ggplot2)
library(ggsci)
library(tidyverse)
library(reshape2)

##############################################################################

gps = read.table('01-samplesGroup.txt', header = T, sep = '\t', stringsAsFactors = F)
rownames(gps) = gps$SampleID
taxon = read.table('01-Taxon.txt',  header = T, sep = '\t', 
                   row.names = 1, stringsAsFactors = F)
otuDF = read.table('01-OTUtable.txt', header = T, sep = '\t',
                   row.names = 1, stringsAsFactors = F)
taxon = taxon[rownames(otuDF), ]

m.otu = matrix(colSums(otuDF), ncol = ncol(otuDF), 
               nrow = nrow(otuDF), byrow = T)
otuDF = otuDF/m.otu

##############################################################################
fams = readLines('01-Top_useFamily.txt')
# fams = c('f_Marinifilaceae', 'f_Marinilabiliaceae','f_Prolixibacteraceae')

otuDF = otuDF[rownames(otuDF)[taxon$Family %in% fams], ]
otuDF = otuDF[rowSums(otuDF) != 0, ]
otuDF$Taxon = taxon[rownames(otuDF), 'Family']
fm.otuDF = aggregate(x = otuDF[, 1:27], by = list(otuDF$Taxon), FUN = sum) 
colnames(fm.otuDF)[1] = 'Taxon'

ml.otuDF = melt(fm.otuDF, 
                measure.vars = colnames(fm.otuDF)[2:28],
                variable.name = 'SampleID')
ml.otuDF$TimeGP = paste0('D', substr(ml.otuDF$SampleID, 2, 3))
ml.otuDF$Taxon = gsub('f_', '', ml.otuDF$Taxon, perl = T)
ml.otuDF$GP = substr(ml.otuDF$SampleID, 1, 1)

p = ggplot(ml.otuDF, aes(x = Taxon, y = value, fill = GP)) + 
  geom_boxplot(size = 0.1) + 
  scale_fill_manual(values = c('#E69F00', 'grey55', '#0072B2') ) + # C, HEPES, Molybdate
  facet_grid(TimeGP~., scales = 'free') +
  labs(y = 'The average relative abundance', x = NULL) + 
  theme(plot.margin = margin(t = 0.1, r = 1, b = 0.1, l = 0.1, "in"),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(), 
        strip.text.y = element_text(size = 10, color = 'black'),
        axis.title =  element_text(size = 10),
        axis.text.y = element_text(size = 8, color = 'black', hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 8, color = 'black', angle = 20, 
                                   face = 'italic', hjust = 1, vjust = 1), 
        legend.position = 'none') 

ggsave('03-pHbufferTop.pdf', p, units = 'in', width = 3, height = 2.6)

##############################################################################

phdf = read.table('01-pHranges.txt', header = T, sep = '\t')

phdf$incuTime = as.numeric(substr(phdf$Time, 2, 3))
# phdf$incuTime = factor(phdf$incuTime,
#                        levels = c(5, 12, 30),
#                        labels = c('D05', 'D12', 'D30'))

head(phdf)

p = ggplot(phdf, aes(x = incuTime, y = pH, group = incuTime)) + 
  geom_boxplot(size = 0.1, outlier.color = 'grey55') + 
  facet_wrap(.~Group, ncol = 1, scales = 'free_x') + 
  labs(y = 'pH', x = NULL) +
  theme(plot.margin = margin(t = 0.05, r = 0.15, b = 0.05, l = 0, "in"),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(), 
        # strip.text.y = element_text(size = 6, color = 'black'),
        strip.text.x = element_blank(),
        axis.title.y =  element_text(size = 8),
        axis.text.y = element_text(size = 6, color = 'black', hjust = 0, vjust = 0.5),
        axis.text.x = element_text(size = 6, color = 'black', angle = 0, 
                                   hjust = 0.5, vjust = 0.5), 
        legend.position = 'none') 

ggsave('03-pHbufferTop_range.pdf', p, units = 'in', width = 1.1, height = 2.5)

##############################################################################


