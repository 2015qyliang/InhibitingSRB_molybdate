# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

#############################################################################

library(microeco)
library(ape)
library(ggplot2)
library(ggsci)
library(reshape2)

#############################################################################

sample_info = read.table('02-samplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('02-OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('02-comTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('02-FastTree.txt')

dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)

dataset$cal_abund()

#############################################################################

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


#############################################################################

t1 = trans_abund$new(dataset = dataset, taxrank = 'Family',
                     ntaxa = 20, groupmean = "Group")


gps = t1$data_abund$Sample
taxs = t1$data_abund$Taxonomy
abund = t1$data_abund$Abundance
newdf = data.frame(abund, taxs, gps)

newdf = newdf[newdf$taxs %in% t1$data_taxanames, ]

# 
newdf.re = dcast(newdf, gps~taxs, value.var = 'abund' )

newdf.re.m = newdf.re[, 2:21]
newdf.re$others = 100 - rowSums(newdf.re.m)
colnames(newdf.re)[1] = 'family'
newdf.re$family = gsub('T', 'D', newdf.re$family)

m.fams = as.data.frame(matrix(colnames(newdf.re), nrow = 1))
colnames(m.fams) = colnames(newdf.re)
m.orders = as.data.frame(matrix(c('family', 1:21), nrow = 1))
colnames(m.orders) = colnames(newdf.re)

fCols = c()
for (cl in c(colorRampPalette(pal_npg()(10))(20), 'grey70')) {
  fCols = append(fCols, paste(col2rgb(cl), sep = '', collapse = ','))
}
m.colrs = as.data.frame(matrix(c('family', fCols), nrow = 1))
colnames(m.colrs) = colnames(newdf.re)
df.head = rbind(m.orders, m.colrs, m.fams)
df.head = cbind(matrix('family', nrow = 3, 2), df.head)

csM = col2rgb(pal_jco()(4))
# cCols = paste(col2rgb('#E69F00'), sep = '', collapse = ',')
# mCols = paste(col2rgb('#0072B2'), sep = '', collapse = ',')
df.orderCols = data.frame(Od = 25:22, 
                          samCols = c(paste(csM[, 1], sep = '', collapse = ','),
                                      paste(csM[, 2], sep = '', collapse = ','),
                                      paste(csM[, 3], sep = '', collapse = ','),
                                      paste(csM[, 4], sep = '', collapse = ',')))
newdf.re = cbind(df.orderCols, newdf.re)

write.table(df.head, '05-CircosDF.txt', append = F, quote = F, sep = '\t',
            row.names = F,col.names = F)

write.table(newdf.re, '05-CircosDF.txt', append = T, quote = F, sep = '\t',
            row.names = F,col.names = F)

# '#E69F00' --C
# '#0072B2' --M
