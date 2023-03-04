
# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##############################################################

library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggsci)
library(grid)

##############################################################
DiffTaxsDF = read.table('01-diffTaxsAmongGroups.txt', header = T, 
                        sep = '\t', stringsAsFactors = F)

cutoff.abund = 0.5 # choose family according to abundance

cgtm.df = DiffTaxsDF[DiffTaxsDF$Log2Change > 0 & DiffTaxsDF$P.adj <= 0.05, ]
cgtm.df = cgtm.df[which(cgtm.df$Cmean >= cutoff.abund), ]
cgtm.df = cgtm.df[!(str_detect(cgtm.df$Taxs, '_norank_') | 
                      str_detect(cgtm.df$Taxs, '_unclassified_')), ]

clem.df = DiffTaxsDF[DiffTaxsDF$Log2Change < 0 & DiffTaxsDF$P.adj <= 0.05, ]
clem.df = clem.df[which(clem.df$Mmean >= cutoff.abund), ]
clem.df = clem.df[!(str_detect(clem.df$Taxs, '_norank_') | 
                      str_detect(clem.df$Taxs, '_unclassified_')), ]

cgtm.df$Taxs = gsub('f_', '', cgtm.df$Taxs)
clem.df$Taxs = gsub('f_', '', clem.df$Taxs)

tax.fac = rev(c( sort(setdiff(unique(cgtm.df$Taxs), unique(clem.df$Taxs))),
                 intersect(unique(cgtm.df$Taxs), unique(clem.df$Taxs)),
                 sort(setdiff(unique(clem.df$Taxs), unique(cgtm.df$Taxs)))))

##############################################################
##############################################################
##############################################################

DiffTaxsDF = read.table('01-diffTaxsAmongGroups.txt', header = T, 
                        sep = '\t', stringsAsFactors = F)

cutoff.abund = 0 # choose family according to abundance

cgtm.df = DiffTaxsDF[DiffTaxsDF$Log2Change > 0 & DiffTaxsDF$P.adj <= 0.05, ]
cgtm.df = cgtm.df[which(cgtm.df$Cmean >= cutoff.abund), ]
cgtm.df = cgtm.df[!(str_detect(cgtm.df$Taxs, '_norank_') | 
                      str_detect(cgtm.df$Taxs, '_unclassified_')), ]

clem.df = DiffTaxsDF[DiffTaxsDF$Log2Change < 0 & DiffTaxsDF$P.adj <= 0.05, ]
clem.df = clem.df[which(clem.df$Mmean >= cutoff.abund), ]
clem.df = clem.df[!(str_detect(clem.df$Taxs, '_norank_') | 
                      str_detect(clem.df$Taxs, '_unclassified_')), ]

cgtm.df$Taxs = gsub('f_', '', cgtm.df$Taxs)
clem.df$Taxs = gsub('f_', '', clem.df$Taxs)

newdf = rbind(cgtm.df, clem.df)

newdf = newdf[newdf$Taxs %in% tax.fac, ]
newdf$Taxs = factor(newdf$Taxs, levels = tax.fac)

newdf$Groups = substr(newdf$Groups, 1, 3)
newdf$Groups = factor(newdf$Groups, levels = sort(unique(newdf$Groups)))
newdf$signif = rep(' ', nrow(newdf))
newdf$signif[which(newdf$P.adj <= 0.001)] = '***'
newdf$signif[which(newdf$P.adj > 0.001 & newdf$P.adj <= 0.01)] = '**'
newdf$signif[which(newdf$P.adj > 0.01 & newdf$P.adj <= 0.05)] = '*'

###############################################################

newdf$Groups = factor(newdf$Groups, levels = c('T30', 'T21', 'T12', 'T05'))

# newdf$Groups = factor(newdf$Groups, levels = rev(c('T30', 'T21', 'T12', 'T05')))

p1 = ggplot(newdf, aes(x = Log2Change, y = Groups, fill = Groups)) + 
  geom_col(position = 'dodge') + 
  geom_text(label = newdf$signif, size = 2.3, nudge_y = -0.45) +
  scale_fill_manual(values = rev(pal_jco()(4))) + 
  facet_grid(Taxs~., scales = 'free_y', space = 'free_y') + 
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6), 
        # axis.text.y = element_text(size = 6), 
        axis.text.y = element_blank(), 
        strip.text.y = element_text(size = 8, angle = 360, hjust = 1, face = 'italic'), 
        strip.background.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key =  element_blank(),
        legend.key.height = unit(1, 'mm'),
        legend.text = element_text(size = 6))

p1

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

newdf = rbind(cgtm.df, clem.df)
newdf = newdf[newdf$Taxs %in% tax.fac, ]

DFmelt = c()
for (gp in unique(newdf$Groups)) {
  gpc = unlist(strsplit(gp, split = '_vs_'))[1]
  gpm = unlist(strsplit(gp, split = '_vs_'))[2]
  gDF = newdf[which(newdf$Groups == gp), ]
  gDF = melt(gDF, measure.vars = c('Cmean', 'Mmean'), 
             variable.name = 'Grp', value.name = 'Abundance')
  gDF$Grp = as.character(gDF$Grp)
  gDF$Grp[which(gDF$Grp == 'Cmean')] = gpc
  gDF$Grp[which(gDF$Grp == 'Mmean')] = gpm
  DFmelt = rbind(DFmelt, gDF)
}

DFmelt$Taxs = factor(DFmelt$Taxs, levels = tax.fac)

DFmelt$TimeGP = substr(DFmelt$Grp, 1, 3)
DFmelt$Grp = substr(DFmelt$Grp, 4, 4)

# DFmelt$Abundance[DFmelt$Abundance < 0.5] = NA

DFmelt$TimeGP = factor(DFmelt$TimeGP, levels = c('T30', 'T21', 'T12', 'T05'))


p2 = ggplot(DFmelt, aes(x = Abundance, y = TimeGP, fill = TimeGP)) + 
  geom_col(position = 'dodge') + 
  scale_fill_manual(values = rev(pal_jco()(4))) + 
  facet_grid(Taxs~Grp, scales = 'free_y', space = 'free_y') + 
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_blank(), 
        strip.text.y = element_text(size = 6, angle = 360),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key =  element_blank(),
        legend.key.height = unit(1, 'mm'),
        legend.text = element_text(size = 6))

p2

###############################################################
# set the layout of result

# pdf(width = 5, height = 4, file = '12-DiffTaxsLog2Abundance.pdf' )
pdf(width = 5, height = 8, file = '03-DiffTaxsLog2Abundance.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1, 9))) 
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1:4))
print(p2, vp = vplayout(1,5:9))
dev.off()

