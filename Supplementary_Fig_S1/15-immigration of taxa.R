# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

############################################

library(ggplot2)
library(ggsci)
library(grid)
library(tidyverse)
library(reshape2)

############################################

comm = read.table('16-OTUtableC.txt', header = T, sep = '\t', 
                   row.names = 1, stringsAsFactors = F)
comm.m = as.matrix(comm)
comm.m[comm.m > 0] = 1
comm.m = comm.m[, sort(colnames(comm.m))]

cutvalue = 5

comm.m.g = comm.m[, 1:20]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.0 = rownames(comm.m.g)

comm.m.g = comm.m[, 21:40]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.05 = rownames(comm.m.g)

comm.m.g = comm.m[, 41:60]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.12 = rownames(comm.m.g)

comm.m.g = comm.m[, 61:80]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.21 = rownames(comm.m.g)

comm.m.g = comm.m[, 81:100]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.30 = rownames(comm.m.g)

# length(taxa.0)
# length(taxa.05)
# length(taxa.12)
# length(taxa.21)
# length(taxa.30)

obs.1 = length(taxa.0)
obs.2 = length(union(taxa.0, taxa.05))
obs.3 = length(union(union(taxa.0, taxa.05), taxa.12))
obs.4 = length(union(union(union(taxa.0, taxa.05), taxa.12), taxa.21))
obs.5 = length(union(union(union(union(taxa.0, taxa.05), taxa.12), taxa.21), taxa.30))

pre.1 = length(taxa.0)
pre.2 = length(setdiff(taxa.05, taxa.0))
pre.3 = length(setdiff(taxa.12, union(taxa.0, taxa.05) ))
pre.4 = length(setdiff(taxa.21, union(union(taxa.0, taxa.05), taxa.12) ))
pre.5 = length(setdiff(taxa.30, union(union(union(taxa.0, taxa.05), taxa.12), taxa.21)))

gp = rep('C', 5)
t.line = c('D00', 'D05', 'D12', 'D21', 'D30')
pre.unobs = c(pre.1, pre.2, pre.3, pre.4, pre.5)
obs = c(obs.1, obs.2, obs.3, obs.4, obs.5)
frac = round(pre.unobs/obs, digits = 2)
cDF = data.frame(gp, t.line, frac, pre.unobs, obs)

############################################

comm = read.table('16-OTUtableM.txt', header = T, sep = '\t', 
                  row.names = 1, stringsAsFactors = F)
comm.m = as.matrix(comm)
comm.m[comm.m > 0] = 1
comm.m = comm.m[, sort(colnames(comm.m))]

comm.m.g = comm.m[, 1:20]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.0 = rownames(comm.m.g)

comm.m.g = comm.m[, 21:40]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.05 = rownames(comm.m.g)

comm.m.g = comm.m[, 41:60]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.12 = rownames(comm.m.g)

comm.m.g = comm.m[, 61:80]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.21 = rownames(comm.m.g)

comm.m.g = comm.m[, 81:100]
comm.m.g = comm.m.g[rowSums(comm.m.g) >= cutvalue, ]
taxa.30 = rownames(comm.m.g)

# length(taxa.0)
# length(taxa.05)
# length(taxa.12)
# length(taxa.21)
# length(taxa.30)

obs.1 = length(taxa.0)
obs.2 = length(union(taxa.0, taxa.05))
obs.3 = length(union(union(taxa.0, taxa.05), taxa.12))
obs.4 = length(union(union(union(taxa.0, taxa.05), taxa.12), taxa.21))
obs.5 = length(union(union(union(union(taxa.0, taxa.05), taxa.12), taxa.21), taxa.30))

pre.1 = length(taxa.0)
pre.2 = length(setdiff(taxa.05, taxa.0))
pre.3 = length(setdiff(taxa.12, union(taxa.0, taxa.05) ))
pre.4 = length(setdiff(taxa.21, union(union(taxa.0, taxa.05), taxa.12) ))
pre.5 = length(setdiff(taxa.30, union(union(union(taxa.0, taxa.05), taxa.12), taxa.21)))

gp = rep('M', 5)
t.line = c('D00', 'D05', 'D12', 'D21', 'D30')
pre.unobs = c(pre.1, pre.2, pre.3, pre.4, pre.5)
obs = c(obs.1, obs.2, obs.3, obs.4, obs.5)
frac = round(pre.unobs/obs, digits = 2)
mDF = data.frame(gp, t.line, frac, pre.unobs, obs)

############################################

DFall = rbind(cDF, mDF)

DFall.melt = melt(DFall, measure.vars = c('frac', 'pre.unobs', 'obs'), 
                  variable.name = 'Index')
DFall.melt$Index = factor(DFall.melt$Index, 
                          levels = c('frac', 'pre.unobs', 'obs'))

############################################

DFall.melt.filt = DFall.melt[DFall.melt$t.line != 'D00' & 
                                 DFall.melt$Index != 'obs', ]

p1 = ggplot(DFall.melt.filt, aes(x = t.line, y = value, group = gp)) + 
  geom_point(aes(color = gp), shape = 16, size = 1, alpha = 0.8) +
  geom_line(aes(color = gp), linetype = 1, size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c('#E69F00', '#0072B2')) +
  facet_wrap(.~Index, nrow = 1, scales = "free_y") +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_line(size = 0.1, linetype = 2, color = 'grey80'),
        strip.background = element_blank(),
        # legend.position = c(0.2, 0.6), 
        legend.position = 'bottom', 
        legend.background = element_blank(),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'))

DFall.melt.filt = DFall.melt[ DFall.melt$Index == 'obs', ]
p2 = ggplot(DFall.melt.filt, aes(x = t.line, y = value, group = gp)) + 
  geom_point(aes(color = gp), shape = 16, size = 1, alpha = 0.8) +
  geom_line(aes(color = gp), linetype = 1, size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c('#E69F00', '#0072B2')) +
  facet_wrap(.~Index, nrow = 1) +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_line(size = 0.1, linetype = 2, color = 'grey80'),
        strip.background = element_blank(),
        # legend.position = c(0.2, 0.6), 
        legend.position = 'bottom', 
        legend.background = element_blank(),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'))

######################################################################

# set the layout of result
pdf(width = 6, height = 1.75, file = paste0('17-ImmigrationTaxa_filt-', cutvalue,'.pdf') )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1, 3))) 
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(1, 3))
dev.off()


