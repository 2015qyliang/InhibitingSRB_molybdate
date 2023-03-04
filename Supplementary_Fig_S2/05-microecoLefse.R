# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##################################################################

library(tidyverse)
library(ggplot2)
library(ggsci)

##################################################################

lefseDF = read.table('04-res_lefse.txt', header = T, sep = '\t', stringsAsFactors = F)
lefseDF = lefseDF[lefseDF$LDA >= 4, ]

colnames(lefseDF)

p = ggplot(lefseDF, aes(x = Taxa, y = LDA, fill = GP)) + 
  geom_col(alpha = 0.7, width = 0.7) + 
  geom_hline(yintercept = 4, colour="grey70", linetype="dashed")+
  scale_fill_manual(values = c('#E69F00', '#0072B2')) + # C, M
  # scale_fill_nejm() + 
  facet_grid(GP~Group, scales = 'free', space = "free") + 
  labs(fill = NULL, x = NULL, y = 'LDA Score') +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8, color="black", hjust = 1, vjust = 0.5, face = 'italic'),
        axis.text.x = element_text(size = 8, color="black", hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none") + 
  coord_flip()

ggsave('06-LefSePlot.pdf', p, units = 'in', width = 6, height = 2.5)


