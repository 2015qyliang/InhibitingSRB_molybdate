# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

#############################################################

library(vegan)
library(ggsci)
library(ggplot2)
library(tidyverse)

#############################################################

diyscale = function(x) {
  x= (x - min(x))/(max(x) - min(x))
  return(x)
}

######################################

otuSampleDF = t(read.table('01-IndexPearson.txt', header = T,sep = '\t',row.names = 1))

df.scale = apply(otuSampleDF, 2, diyscale)

nmds = metaMDS(vegdist(df.scale, method = 'bray', binary = F), 
               k = 2, trymax = 100, wascores = TRUE)
nmds$stress

data.scores = as.data.frame(scores(nmds))
dt.scores = data.frame(data.scores, Group = rownames(otuSampleDF))

dt.scores.m = dt.scores[dt.scores$Group == 'D00C', ]
dt.scores.m$Group = 'D00M'

newdf = rbind(dt.scores, dt.scores.m)

newdf$GP = substr(newdf$Group, 4, 4)

######################################
######################################

c.newdf = newdf[newdf$GP == 'C', ]
m.newdf = newdf[newdf$GP == 'M', ]
m.newdf = m.newdf[c('D00C1', 'D05M', 'D12M', 'D21M', 'D30M'), ]

c.newdf$NMDS1x = c(c.newdf$NMDS1[2:length(c.newdf$NMDS1)], c.newdf$NMDS1[length(c.newdf$NMDS1)])
c.newdf$NMDS2y = c(c.newdf$NMDS2[2:length(c.newdf$NMDS2)], NA)

m.newdf$NMDS1x = c(m.newdf$NMDS1[2:length(m.newdf$NMDS1)], m.newdf$NMDS1[length(m.newdf$NMDS1)])
m.newdf$NMDS2y = c(m.newdf$NMDS2[2:length(m.newdf$NMDS2)], NA)


nmds.p = ggplot(newdf, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = c.newdf, aes(x = NMDS1, y = NMDS2),
             color = '#E69F00', size = 3, alpha = 0.7, shape = 16) +
  geom_text(data = c.newdf, aes(x = NMDS1, y = NMDS2, label = Group), 
            hjust = 0.5, nudge_x = 0, vjust = 0, nudge_y = 0, size = 2) + 
  geom_segment(data = c.newdf, aes(x = NMDS1, y = NMDS2, xend = NMDS1x , yend = NMDS2y), 
               arrow = arrow(ends = 'last', angle = 20, length = unit(0.1, "in")), 
               colour = '#E69F00', size = .5, alpha = 0.3) + 
  geom_point(data = m.newdf, aes(x = NMDS1, y = NMDS2),
             color = '#0072B2', size = 3, alpha = 0.7, shape = 17) +
  geom_text(data = m.newdf, aes(x = NMDS1, y = NMDS2, label = Group), 
            hjust = 0.5, nudge_x = 0, vjust = 0, nudge_y = 0, size = 2) + 
  geom_segment(data = m.newdf, aes(x = NMDS1, y = NMDS2, xend = NMDS1x , yend = NMDS2y), 
               arrow = arrow(ends = 'last', angle = 20, length = unit(0.1, "in")), 
               colour = '#0072B2', size = .5, alpha = 0.3) + 
  labs(caption = paste0('Stress: ', nmds$stress)) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
  geom_vline(xintercept=0, colour="grey60", linetype="dashed")+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 8), 
        plot.caption =  element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 8),
        legend.position = "none")  +
  guides(color = guide_legend(nrow = 1,
                              keywidth = unit(1, 'mm'),
                              keyheight = unit(1, 'mm'),
                              label.position = 'right',
                              title.position = 'top',
                              direction = 'horizonal'))
ggsave('03-NMDS-networkParameters_Pearson.pdf', 
       nmds.p, units = 'in', width = 3, height = 3)

#############################################################
#############################################################
#############################################################
#############################################################


otuSampleDF = t(read.table('01-IndexSpearman.txt', header = T,sep = '\t',row.names = 1))

df.scale = apply(otuSampleDF, 2, diyscale)

nmds = metaMDS(vegdist(df.scale, method = 'bray', binary = F), 
               k = 2, trymax = 100, wascores = TRUE)
nmds$stress

data.scores = as.data.frame(scores(nmds))
dt.scores = data.frame(data.scores, Group = rownames(otuSampleDF))

dt.scores.m = dt.scores[dt.scores$Group == 'D00C', ]
dt.scores.m$Group = 'D00M'

newdf = rbind(dt.scores, dt.scores.m)

newdf$GP = substr(newdf$Group, 4, 4)

######################################
######################################

c.newdf = newdf[newdf$GP == 'C', ]
m.newdf = newdf[newdf$GP == 'M', ]
m.newdf = m.newdf[c('D00C1', 'D05M', 'D12M', 'D21M', 'D30M'), ]

c.newdf$NMDS1x = c(c.newdf$NMDS1[2:length(c.newdf$NMDS1)], c.newdf$NMDS1[length(c.newdf$NMDS1)])
c.newdf$NMDS2y = c(c.newdf$NMDS2[2:length(c.newdf$NMDS2)], NA)

m.newdf$NMDS1x = c(m.newdf$NMDS1[2:length(m.newdf$NMDS1)], m.newdf$NMDS1[length(m.newdf$NMDS1)])
m.newdf$NMDS2y = c(m.newdf$NMDS2[2:length(m.newdf$NMDS2)], NA)


nmds.p = ggplot(newdf, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = c.newdf, aes(x = NMDS1, y = NMDS2),
             color = '#E69F00', size = 3, alpha = 0.7, shape = 16) +
  geom_text(data = c.newdf, aes(x = NMDS1, y = NMDS2, label = Group), 
            hjust = 0.5, nudge_x = 0, vjust = 0, nudge_y = 0, size = 2) + 
  geom_segment(data = c.newdf, aes(x = NMDS1, y = NMDS2, xend = NMDS1x , yend = NMDS2y), 
               arrow = arrow(ends = 'last', angle = 20, length = unit(0.1, "in")), 
               colour = '#E69F00', size = .5, alpha = 0.3) + 
  geom_point(data = m.newdf, aes(x = NMDS1, y = NMDS2),
             color = '#0072B2', size = 3, alpha = 0.7, shape = 17) +
  geom_text(data = m.newdf, aes(x = NMDS1, y = NMDS2, label = Group), 
            hjust = 0.5, nudge_x = 0, vjust = 0, nudge_y = 0, size = 2) + 
  geom_segment(data = m.newdf, aes(x = NMDS1, y = NMDS2, xend = NMDS1x , yend = NMDS2y), 
               arrow = arrow(ends = 'last', angle = 20, length = unit(0.1, "in")), 
               colour = '#0072B2', size = .5, alpha = 0.3) + 
  labs(caption = paste0('Stress: ', nmds$stress)) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
  geom_vline(xintercept=0, colour="grey60", linetype="dashed")+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 8), 
        plot.caption =  element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 8),
        legend.position = "none")  +
  guides(color = guide_legend(nrow = 1,
                              keywidth = unit(1, 'mm'),
                              keyheight = unit(1, 'mm'),
                              label.position = 'right',
                              title.position = 'top',
                              direction = 'horizonal'))
ggsave('03-NMDS-networkParameters_Spearman.pdf', 
       nmds.p, units = 'in', width = 3, height = 3)

