# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##################################################################

library(ggplot2)
library(grid)
library(doBy)

# c('#E69F00', '#0072B2')

##########################################

pHdf = read.table('01-pH_Optimum.txt', header = T, sep = '\t', stringsAsFactors = F)

p1 = ggplot(pHdf, aes(x = family, y = optimum)) + 
  geom_boxplot(size = 0.25, shape = 16, outlier.size = 1, outlier.color = 'grey75') +
  geom_jitter(alpha = 0.4, color = 'grey50', shape = 16, width = 0.25) + 
  geom_hline(yintercept = 7.426, color = '#E69F00', size = 0.3) +
  geom_rect(aes(ymin= 7.21, ymax= 7.64, xmin=0, xmax=5),
            fill='#E69F00', alpha = 0.01) +
  geom_hline(yintercept = 6.435, color = '#0072B2', size = 0.3) +
  geom_rect(aes(ymin= 6.2, ymax= 6.67, xmin=0, xmax=5),
            fill='#0072B2', alpha = 0.01) +
  scale_y_continuous(limits =  c(6, 10)) + 
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black',  
                                   hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = 'black', angle = 90,
                                   hjust = 0.5, vjust = 0.5))


##########################################

pHdf = read.table('01-pH_Minimum.txt', header = T, sep = '\t', stringsAsFactors = F)

p2 = ggplot(pHdf, aes(x = family, y = minimum)) + 
  geom_boxplot(size = 0.25, shape = 16, outlier.size = 1, outlier.color = 'grey75') +
  geom_jitter(alpha = 0.4, color = 'grey50', shape = 16, width = 0.25) + 
  geom_hline(yintercept = 6.435, color = '#0072B2', size = 0.3) +
  geom_rect(aes(ymin= 6.2, ymax= 6.67, xmin=0, xmax=5),
            fill='#0072B2', alpha = 0.01) +
  geom_hline(yintercept = 7.426, color = '#E69F00', size = 0.3) +
  geom_rect(aes(ymin= 7.21, ymax= 7.64, xmin=0, xmax=5),
            fill='#E69F00', alpha = 0.01) +
  scale_y_continuous(limits =  c(4, 9) ) + 
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black',  hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = 'black', angle = 90,
                                   hjust = 0.5, vjust = 0.5))

##########################################

# set the layout of result
pdf(width = 4.6, height = 4, file = '03-pHrange.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(2, 1))) 
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))
dev.off()

####################################################################################




