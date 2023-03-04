# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##############################################################################

library(ggplot2)
library(ggsci)
library(tidyverse)
library(reshape2)

########################################################

fns0 = list.files('OD600/', 'T0')
fns28 = list.files('OD600/', 'T54')


########################################################

growthDF = data.frame()
for (fi in 1:13) {
  m0 = as.matrix(read.table(file.path('OD600/', fns0[fi]), header = F, sep = '\t'))
  m1 = as.matrix(read.table(file.path('OD600/', fns28[fi]), header = F, sep = '\t'))
  m.growth = as.data.frame(m1 - m0)
  m.growth = m.growth[1:8, ]
  m.growth$media = c('1CK', '2Pro', '3Ori', '4Met', '5B7', '6B12', '7Cys', '8Iso')
  m.growth$Species = rep(gsub('.txt', '', str_split(fns0[fi], '-', simplify = T)[1, 2]), 4)
  growthDF = rbind(growthDF, 
                   melt(m.growth, measure.names = colnames(m.growth)[1:5], value.name = 'Growth'))
}

########################################################

colnames(growthDF)
head(growthDF)

media = Species = Growth = Gsd = c()
for (sp in as.character(unique(growthDF$Species))) {
  for (md in unique(growthDF$media)) {
    spmdDF = growthDF[growthDF$Species == sp & growthDF$media == md, ]
    out.test1 = boxplot.stats(spmdDF$Growth)
    if (length(out.test1$out) != 0) {
      spmdDF.f = spmdDF[spmdDF$Growth != out.test1$out, ]
      media = append(media, md)
      Species = append(Species, sp)
      Growth = append(Growth, mean(spmdDF.f$Growth))
      Gsd = append(Gsd, sd(spmdDF.f$Growth))
    } else {
      spmdDF.f = spmdDF
      media = append(media, md)
      Species = append(Species, sp)
      Growth = append(Growth, mean(spmdDF.f$Growth))
      Gsd = append(Gsd, sd(spmdDF.f$Growth))
    }
  }
}

newdf = data.frame(media, Species, Growth, Gsd)
newdf$Species = factor(newdf$Species,
                       levels = c('S1','S2','S3','S4','S5','S6',
                                  'S7','S8','S9','S10','S11','S12','S13'))

########################################################

p = ggplot(newdf, aes(x = Species, y = Growth, fill = media)) + 
  geom_col(position = 'dodge') + 
  geom_errorbar(aes(ymin = Growth-Gsd, ymax = Growth+Gsd), size = 0.2, 
                width = 0.5, position = position_dodge(0.9)) + 
  scale_fill_npg() +  
  labs(y = 'OD600', x = NULL) + theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, color = 'black'), 
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5), 
        legend.key.height = unit(0.05, 'in'), 
        legend.key.width = unit(0.1, 'in'),
        legend.position = 'bottom')
ggsave('05-OD600_growth.pdf', p, units = 'in', width = 6, height = 3)

########################################################
########################################################
########################################################

#############
# anova TEST

# mark Star: TtestS < 0, TtestP < 0.05

spvec = c('S1','S2','S3','S4','S5','S6',
          'S7','S8','S9','S10','S11','S12','S13')
pr.comb = combn(c('1CK', '2Pro', '3Ori', '4Met', '5B7', '6B12', '7Cys', '8Iso'), 2)

write('Species\tM1\tM2\tTtestS\tTtestP', append = F, '06-Ttest.txt' )

for (sp in spvec) {
  spDF = growthDF[growthDF$Species == sp, ]
  # write('Species\M1\tM2\tTtestS\tTtestP', append = F, 
  #       paste0('12-ANOVA_TurkeyHSDtest_', sp, '.txt') )
  for (pr in 1:ncol(pr.comb)) {
    prcs = pr.comb[, pr]
    spDF.f1 = spDF[spDF$media == prcs[1], ]
    spDF.f2 = spDF[spDF$media == prcs[2], ]
    out.test1 = boxplot.stats(spDF.f1$Growth)
    out.test2 = boxplot.stats(spDF.f2$Growth)
    if (length(out.test1$out) != 0) {
      spDF.f1 = spDF.f1[spDF.f1$Growth != out.test1$out, ]
    }
    if (length(out.test2$out) != 0) {
      spDF.f2 = spDF.f2[spDF.f2$Growth != out.test2$out, ]
    }
    spDF.f = rbind(spDF.f1, spDF.f2)
    res.t = t.test(Growth ~ media, data = spDF.f)
    ln.res = paste0(sp, '\t', paste0(prcs, collapse = '\t'), '\t',
                    res.t$statistic, '\t', res.t$p.value)
    # write(ln.res, append = T, paste0('12-Ttest_', sp, '.txt') )
    write(ln.res, append = T, '06-Ttest.txt' )
  }
}

