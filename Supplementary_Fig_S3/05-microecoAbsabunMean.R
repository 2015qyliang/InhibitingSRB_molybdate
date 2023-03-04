# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

#############################################################################

library(microeco)
library(ape)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

#############################################################################

sample_info = read.table('04-samplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('04-OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('04-comTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('04-FastTree.txt')
dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)

dataset$cal_abund()

# #############################################################################
diyplotBar = function(dataset, tax, numRank) {
  t1 = trans_abund$new(dataset = dataset, taxrank = tax,
                       ntaxa = numRank, groupmean = "Group")
  p = t1$plot_bar(color_values  = colorRampPalette(pal_npg()(10))(numRank),
                  others_color  = "grey70",
                  xtext_keep = FALSE, legend_text_italic = T) + 
    theme(panel.background = element_rect(fill="white",color="black"),
          panel.grid = element_blank(), 
          axis.title.y.left =  element_text(size = 10),
          axis.text.x = element_text(size = 8, color = 'black', 
                                     hjust = 0.5, vjust = 0.5), 
          axis.text.y = element_text(size = 8, color = 'black'), 
          legend.key = element_blank(), 
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) + 
    guides(fill = guide_legend( ncol = 3,
                                keywidth = unit(1, 'mm'),
                                keyheight = unit(1, 'mm'),
                                label.hjust = 0,
                                label.vjust = 0.5,
                                label.position = 'right',
                                title.position = 'top', 
                                direction = 'horizonal'))
  return(p)
}

# #############################################################################

newtax = dataset$taxa_abund
for (i in 2:length(newtax)) {
  df1 = newtax[[i]]
  newRN = c()
  for (rn in rownames(df1)) {
    rn = strsplit(rn, split =  '|', fixed = T)[[1]][i]
    newRN = append(newRN, rn)
  }
  rownames(newtax[[i]]) = newRN
}
dataset$taxa_abund = newtax

# #############################################################################

tax = 'Family'
dpb = diyplotBar(dataset, tax, 20)
ggsave(paste0('06-bar-', tax,'-Legend.pdf'),
       dpb, units = 'in', width = 6, height = 3.5)


t1 = trans_abund$new(dataset = dataset, taxrank = tax,
                     ntaxa = 20, groupmean = "Group")

write(t1$use_taxanames, '07-Top20_useFamily.txt')
