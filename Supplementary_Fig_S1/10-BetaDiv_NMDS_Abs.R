# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(vegan)
library(ggsci)
library(ggplot2)
library(tidyverse)

#############################################################

dclassDF = read.table('02-samplesGroup.txt', header = T, 
                      sep = '\t', stringsAsFactors = F)

#############################################################

otuSampleDF = t(read.table('11-AbsTableFamily.txt', header = T,sep = '\t',row.names = 1))
nmds = metaMDS(vegdist(otuSampleDF, method = 'jaccard', binary = T), 
               k = 2, trymax = 100, wascores = TRUE)
info.stress = round(nmds$stress, digits = 2)
data.scores = as.data.frame(scores(nmds))
newdf = data.frame(data.scores, dclass = dclassDF$Group, 
                   TaxLevel = rep('Family', nrow(data.scores)))

write.table(newdf, '12-NMDS_Family.txt', append = F, quote = F,
            sep = '\t', row.names = F, col.names = T)

#############################################################

c.newdf = newdf[which(newdf$dclass %in% c('T00C', 'T05C', 'T12C', 'T21C', 'T30C')), ]
c.newdf$dclass = factor(c.newdf$dclass,
                        labels = c('D00C', 'D05C', 'D12C', 'D21C', 'D30C'),
                        levels = c('T00C', 'T05C', 'T12C', 'T21C', 'T30C'))
m.newdf = newdf[which(newdf$dclass %in% c('T05M', 'T12M', 'T21M', 'T30M')), ]
m.newdf$dclass = factor(m.newdf$dclass,
                        labels = c('D05M', 'D12M', 'D21M', 'D30M'),
                        levels =  c('T05M', 'T12M', 'T21M', 'T30M'))

nmds.p = ggplot(newdf, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = c.newdf, aes(x = NMDS1, y = NMDS2, color = dclass),
             size = 1.5, alpha = 0.7, shape = 16) +
  geom_point(data = m.newdf, aes(x = NMDS1, y = NMDS2, color = dclass),
             size = 1.5, alpha = 0.7, shape = 17) +
  labs(caption = paste0('Stress: ', paste(info.stress, collapse = ', '))) +
  scale_colour_manual(values = c(pal_npg()(5), pal_npg()(5)[2:5]),
                      breaks = c('D00C', 'D05C', 'D12C', 'D21C', 'D30C',
                                 'D05M', 'D12M', 'D21M', 'D30M')) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
  geom_vline(xintercept=0, colour="grey60", linetype="dashed")+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 8), 
        plot.caption =  element_text(size = 8), 
        axis.title = element_text(size = 8),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 8, color = 'black'),
        strip.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_blank())  +
  guides(color = guide_legend(nrow = 1,
                              keywidth = unit(1, 'mm'),
                              keyheight = unit(1, 'mm'),
                              label.position = 'right',
                              title.position = 'top',
                              direction = 'horizonal'))
ggsave('13-BetaNMDS_Family.pdf', 
       nmds.p, units = 'in', width = 2, height = 1.8)


#############################################################
#############################################################
#############################################################

compares = list(c('T00C', 'T05C'), c('T00C', 'T12C'), c('T00C', 'T21C'), 
                c('T00C', 'T30C'), c('T00C', 'T05M'), c('T00C', 'T12M'),
                c('T00C', 'T21M'), c('T00C', 'T30M'), c('T05C', 'T05M'),
                c('T12C', 'T12M'), c('T21C', 'T21M'), c('T30C', 'T30M'))

###############################

otuSampleDF = t(read.table('11-AbsTableFamily.txt', header = T,sep = '\t',row.names = 1))
otuSampleDF = as.matrix(otuSampleDF)
otuSampleDF[otuSampleDF != 0] = 1

adonis.F = adonis.p = anosim.r = anosim.p = mrpp.delta = mrpp.p = gpPairs = c()
for (i in 1:length(compares)) {
  choseSamples = str_detect(rownames(otuSampleDF), compares[[i]][1]) | 
    str_detect(rownames(otuSampleDF), compares[[i]][2])
  choseDF = otuSampleDF[choseSamples, ]
  choseDF = choseDF[, colSums(choseDF) != 0]
  res.adonis = adonis(choseDF ~ substr(rownames(choseDF), 1, 4),
                      permutations = 999, method = 'jaccard')
  res.anosim = anosim(choseDF, substr(rownames(choseDF), 1, 4), 
                      permutations = 999, distance = 'jaccard')
  res.mrpp = mrpp(choseDF, substr(rownames(choseDF), 1, 4), 
                  permutations = 999, distance = 'jaccard')
  adonis.F = append(adonis.F, res.adonis$aov.tab$F.Model[1])
  adonis.p = append(adonis.p, res.adonis$aov.tab$`Pr(>F)`[1])
  anosim.r = append(anosim.r, res.anosim$statistic)
  anosim.p = append(anosim.p, res.anosim$signif)
  mrpp.delta = append(mrpp.delta, res.mrpp$delta)
  mrpp.p = append(mrpp.p, res.mrpp$Pvalue)
  gpPairs = append(gpPairs, paste(compares[[i]], collapse = '_vs_'))
}

write.table(data.frame(gpPairs, adonis.F, adonis.p, anosim.r, anosim.p, mrpp.delta, mrpp.p), 
            '14-CompositionTest_Family.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)


