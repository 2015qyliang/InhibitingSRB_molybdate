# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

########################################################################

library(tidyverse)
library(ggplot2)
library(ggsci)

########################################################################

decoderDF = read.table('01-DecoderFresh.txt', header = T, sep = '\t', stringsAsFactors = F)
decoderDF$sulfur.assimilation[decoderDF$sulfur.assimilation != 1] = 0

vec.fams = c()
vec.path = c()
vec.per = c()
for (fam in unique(decoderDF$Family)) {
  fam.all = as.integer(table(decoderDF$Family)[fam])
  for (pathw in colnames(decoderDF)[2:6]) {
    ct = sum(decoderDF[, pathw] == 1 & decoderDF$Family == fam)
    vec.fams = append(vec.fams, fam)
    vec.path = append(vec.path, pathw)
    vec.per = append(vec.per, round(ct/fam.all, digits = 3))
  }
}

tDF = data.frame(table(decoderDF$Family))
t.labs = paste0(tDF$Var1, ' (', tDF$Freq, ')')
names(t.labs) = tDF$Var1

newDF = data.frame(family = vec.fams, pathway = vec.path, Percent = vec.per*100)
# newDF = newDF[newDF$family %in% sort(unique(DiffTaxsDF$Family), F), ]

newDF$family = factor(newDF$family, labels = t.labs[sort(unique(newDF$family))], 
                      levels = sort(unique(newDF$family)))

# newDF$family = factor(newDF$family, labels = t.labs[sort(unique(newDF$family))], 
#                       levels = sort(c(unique(newDF$family), "f_Nitrincolaceae")))

########################################################################

p = ggplot(newDF, aes(x = family, y = Percent, fill = pathway)) + 
  geom_col(position = 'dodge') + 
  scale_fill_nejm(alpha = 0.85) + 
  labs(x = NULL, y = 'The percentage of\n completed pathway (%)') + 
  theme_bw()+
  theme(plot.margin = unit(c(0.1, 0.1,0, 0.35),'in'),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 10, colour = 'black'), 
        axis.text.x = element_text(size = 8, angle = 45, 
                                   colour = 'black', face = 'italic', 
                                   hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, angle = 90), 
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key =  element_blank(),
        legend.key.width = unit(2, 'mm'),
        legend.key.height = unit(2, 'mm'),
        legend.text = element_text(size = 8))

ggsave('03-DecoderKO.pdf', p, width = 6, height = 3 )

