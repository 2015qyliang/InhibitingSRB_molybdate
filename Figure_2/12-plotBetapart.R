# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(ggplot2)
library(ggsci)
library(grid)
library(viridis)
library(tidyverse)
library(reshape2)
library(vegan)


################################################

distmelt = function(d){
  d.melt = d %>% as.matrix %>% as.data.frame %>% tibble::rownames_to_column() %>% melt()
  d.melt$variable = as.character(d.melt$variable)
  return(d.melt[d.melt$variable > d.melt$rowname, ])
}

# Better combined graphs
plot_dm = function(df, ylab){
  return(
    ggplot(df, aes(y = variable, x = rowname, fill = value)) + geom_raster() +
      theme(panel.background = element_rect(fill="white", color="white"),
            strip.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = c(0.98, 0.4), 
            legend.background = element_blank(), 
            legend.key.width = unit(1.5, 'mm'), 
            legend.key = element_blank(),
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6)) +
      labs(y = ylab, fill = "Binary\nJaccard\nDistance") +
      facet_grid(var_cat ~ Type + row_cat, scales = "free", switch = "y")
  )
}

###########################################################

# Add columns listing the categories that the samples are in. Match them using grep.
fix_rows = function(df){
  df$Type = factor(df$Type, levels = c("Total", "Nestedness", "Turnover"))
  df$row_cat = 'D00'
  df$row_cat[grepl("T05C", df$rowname, fixed = T)] = 'D05'
  df$row_cat[grepl("T12C", df$rowname, fixed = T)] = 'D12'
  df$row_cat[grepl("T21C", df$rowname, fixed = T)] = 'D21'
  df$row_cat[grepl("T30C", df$rowname, fixed = T)] = 'D30'
  df$row_cat = factor(df$row_cat, levels = c('D00', 'D05', 'D12','D21','D30'))
  df$var_cat = 'D00'
  df$var_cat[grepl("T05C", df$variable, fixed = T)] = 'D05'
  df$var_cat[grepl("T12C", df$variable, fixed = T)] = 'D12'
  df$var_cat[grepl("T21C", df$variable, fixed = T)] = 'D21'
  df$var_cat[grepl("T30C", df$variable, fixed = T)] = 'D30'
  df$var_cat = factor(df$var_cat, levels = rev(c('D00', 'D05', 'D12','D21','D30')))
  return(df)
}

bp = readRDS('11-bpC.rds')
jaccard.melt = bp$beta.jac %>% distmelt()
turnover.melt = bp$beta.jtu %>% distmelt()
nestedness.melt = bp$beta.jne %>% distmelt()

# Merge for combined graph
jaccard.melt$Type = "Total"
turnover.melt$Type = "Turnover"
nestedness.melt$Type = "Nestedness"
all.melt = rbind(jaccard.melt, turnover.melt, nestedness.melt)

all.melt = all.melt %>% fix_rows

p1 = all.melt %>% plot_dm(ylab = "Enrichment Days (Control)") + scale_fill_viridis(limits=c(0,1)) 

ggsave('13-betapart_C.pdf', p1, width = 6, height = 3)

###########################################################

# Add columns listing the categories that the samples are in. Match them using grep.
fix_rows = function(df){
  df$Type = factor(df$Type, levels = c("Total", "Nestedness", "Turnover"))
  df$row_cat = 'D00'
  df$row_cat[grepl("T05M", df$rowname, fixed = T)] = 'D05'
  df$row_cat[grepl("T12M", df$rowname, fixed = T)] = 'D12'
  df$row_cat[grepl("T21M", df$rowname, fixed = T)] = 'D21'
  df$row_cat[grepl("T30M", df$rowname, fixed = T)] = 'D30'
  df$row_cat = factor(df$row_cat, levels = c('D00', 'D05', 'D12','D21','D30'))
  df$var_cat = 'D00'
  df$var_cat[grepl("T05M", df$variable, fixed = T)] = 'D05'
  df$var_cat[grepl("T12M", df$variable, fixed = T)] = 'D12'
  df$var_cat[grepl("T21M", df$variable, fixed = T)] = 'D21'
  df$var_cat[grepl("T30M", df$variable, fixed = T)] = 'D30'
  df$var_cat = factor(df$var_cat, levels = rev(c('D00', 'D05', 'D12','D21','D30')))
  return(df)
}

bp = readRDS('11-bpM.rds')
jaccard.melt = bp$beta.jac %>% distmelt()
turnover.melt = bp$beta.jtu %>% distmelt()
nestedness.melt = bp$beta.jne %>% distmelt()

# Merge for combined graph
jaccard.melt$Type = "Total"
turnover.melt$Type = "Turnover"
nestedness.melt$Type = "Nestedness"
all.melt = rbind(jaccard.melt, turnover.melt, nestedness.melt)

all.melt = all.melt %>% fix_rows

p2 = all.melt %>% plot_dm(ylab = "Enrichment Days (Treat)") + scale_fill_viridis(limits=c(0,1))

ggsave('13-betapart_M.pdf', p2, width = 6, height = 3)

###########################################################


# set the layout of result
pdf(width = 6.5, height = 5, file = '14-betapart.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(2, 1))) 
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))
dev.off()


###########################################################
# summary Info

bp = readRDS('11-bpC.rds')

groupDF = read.table('11-OTUtableC.txt', header = T, sep = '\t', stringsAsFactors = F)
group = paste0('D', substr(colnames(groupDF)[2:ncol(groupDF)], 2, 3))

aovResult = adonis(bp$beta.jac ~ group)$aov.tab
jac.r = aovResult$R2[1]
jac.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jne ~ group)$aov.tab
nes.r = aovResult$R2[1]
nes.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jtu ~ group)$aov.tab
tur.r = aovResult$R2[1]
tur.p = aovResult$`Pr(>F)`[1]

c.values = c(jac.r, jac.p, nes.r, nes.p, tur.r, tur.p)

#########################

bp = readRDS('11-bpM.rds')

groupDF = read.table('11-OTUtableM.txt', header = T, sep = '\t', stringsAsFactors = F)
group = paste0('D', substr(colnames(groupDF)[2:ncol(groupDF)], 2, 3))

aovResult = adonis(bp$beta.jac ~ group)$aov.tab
jac.r = aovResult$R2[1]
jac.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jne ~ group)$aov.tab
nes.r = aovResult$R2[1]
nes.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jtu ~ group)$aov.tab
tur.r = aovResult$R2[1]
tur.p = aovResult$`Pr(>F)`[1]
m.values = c(jac.r, jac.p, nes.r, nes.p, tur.r, tur.p)

#########################

indexs = c('jac.r', 'jac.p', 'nes.r', 'nes.p', 'tur.r', 'tur.p')
newdf = data.frame(indexs, c.values, m.values)

write.table(newdf, '15-adnois.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

