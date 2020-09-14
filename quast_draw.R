library(ggplot2)
library(reshape2)

qdata = as.data.frame(t(as.matrix(read.table('report.tsv', 
                                             row.names=1, sep='\t',
                                             quote='"'))))
qdata[, 'Complete BUSCOs'] = c(1769, 2125, 1657, 2125, 2126)
qdata[, 'Fragmented BUSCOs'] = c(225, 4, 286, 3, 2)
qdata[, 'Missing BUSCOs'] = c(143, 9, 194, 9, 9)

qual_melted = melt(qdata, id.vars="Assembly")
qual_melted$value = as.numeric(qual_melted$value)

tpl = qual_melted[qual_melted$variable %in% c('N50',
                                              'Num. of contigs (>= 0 bp)',
                                              'Num. of mismatches per 100 kbp',
                                              'Num. of indels per 100 kbp',
                                              'Complete BUSCOs',
                                              'Fragmented BUSCOs',
                                              'Missing BUSCOs'
                                              ),]

ggplot(tpl, aes(x=Assembly, y=value, fill=Assembly)) + 
  geom_bar(stat='identity', col='black') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~variable, scales="free", nrow=1) +
  guides(fill=F)


#busci = read.table('busco.tsv', header=T, sep='\t')
#ggplot(busci, aes(x=Strain, y=Count, fill=Type)) + 
#  geom_bar(stat='identity', col='black') +
#  theme_bw() + coord_flip() + guides(fill=F)
