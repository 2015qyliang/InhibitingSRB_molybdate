
# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

######################################################

library(RADanalysis)
library(ggplot2)
library(grid)
library(reshape2)

rankNum = 50
otuDF = read.table('02-OTUtable.txt', header = T, sep = '\t', 
                   row.names = 1, stringsAsFactors = F)
nrads = RADnormalization_matrix(input = otuDF, max_rank = rankNum,
                                average_over = 20, sample_in_row = F)
nrads.m = as.data.frame(nrads$norm_matrix)
colnames(nrads.m) = 1:rankNum
nrads.m$sampleIDs = colnames(otuDF)

nrads.df = melt(nrads.m, id.vars = 'sampleIDs', variable.name = 'nrad', 
                measure.vars = colnames(nrads.m)[1:rankNum])
nrads.df$TimeGroup = substr(nrads.df$sampleIDs, 1, 3)
nrads.df$TimeGroup = factor(nrads.df$TimeGroup,
                            labels = c('D00', 'D05', 'D12', 'D21', 'D30'), 
                            levels = c('T00', 'T05', 'T12', 'T21', 'T30'))
nrads.df$Group = substr(nrads.df$sampleIDs, 4, 4)
nrads.df$Group = factor(nrads.df$Group, levels = c('C', 'M'), 
                        labels = c('Control', 'Treat'))

p1 = ggplot(nrads.df, aes(x = nrad, y = value, color = Group, group = sampleIDs)) + 
  geom_line(alpha = 0.3, size = 0.2) + 
  scale_color_manual(values = c('#E69F00', '#0072B2')) +
  # geom_smooth(method = 'lm', linetype = 2, size = 0.1, se = F, alpha = 0.1) +
  # facet_grid(Group~TimeGroup) + 
  facet_wrap(.~TimeGroup, nrow = 1) + 
  labs(x = 'Rank', y = 'Relative abundance') + 
  scale_y_log10(limits = c(0.001, 1)) + 
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 6, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 0),
        strip.text = element_text(size = 10, color = 'black'),
        strip.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_blank())

######################################################

slope = c()
for (gp in unique(nrads.df$sampleIDs)) {
  n1 = nrads.df[nrads.df$sampleIDs == gp, ]
  n1$nrad = as.integer(n1$nrad)
  lm.rst = lm(log10(n1$value)~n1$nrad)
  sl = lm.rst$coefficients[2]
  slope = append(slope, sl)
}
slopeDF = data.frame(sampleIDs = unique(nrads.df$sampleIDs), slope)
slopeDF$TimeGroup = substr(slopeDF$sampleIDs, 1, 3)
slopeDF$TimeGroup = factor(slopeDF$TimeGroup,
                            labels = c('D00', 'D05', 'D12', 'D21', 'D30'), 
                            levels = c('T00', 'T05', 'T12', 'T21', 'T30'))
slopeDF$Group = substr(slopeDF$sampleIDs, 4, 4)
# slopeDF$Group = factor(slopeDF$Group, levels = c('C', 'M'), 
#                         labels = c('Control', 'Treat'))
slopeDF$gp = paste0('D', substr(slopeDF$sampleIDs, 2, 4))

set.seed(618)
p2 = ggplot(slopeDF, aes(x = gp, y = slope, color = Group)) + 
  geom_boxplot(alpha = 0.3, size = 0.1) + 
  geom_jitter(width = 0.25, alpha = 0.4) + 
  scale_color_manual(values = c('#E69F00', '#0072B2')) +
  facet_grid(.~TimeGroup, scales = 'free_x', space = 'free') + 
  labs(x = NULL, y = 'Slope') +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(), 
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 6, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 10, color = 'black'),
        strip.background = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_blank())

######################################################
# set the layout of result
pdf(width = 8, height = 4, file = '08-RADs_slope.pdf' )
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
print(p1, vp = vplayout(1:3, 1))
print(p2, vp = vplayout(4:5, 1))
dev.off()

######################################################

cmSlopeDF = slopeDF[slopeDF$gp != 'D00C', ]
head(cmSlopeDF)

for (tgp in unique(cmSlopeDF$TimeGroup)) {
  df1 = cmSlopeDF[cmSlopeDF$TimeGroup == tgp, ]
  posthoc = TukeyHSD(aov(slope ~ gp, data = df1), 'gp', conf.level=0.95)
  write.table(data.frame(Compares = rownames(posthoc$gp), posthoc$gp), 
              paste0('09_ANOVA_TurkeyHSD_', tgp, '.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}

