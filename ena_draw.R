# ENA picture

library(lattice)
library(colorRamps)

ena = as.matrix(read.table('ena.tsv', sep='\t', row.names=1, header=T))

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)

drawCoolHM = function(df){
  e = round(df, digits=2)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(-0, 70, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

pdf('/media/barbitoff/DATA/Working issues/FizGen/Translation and Genome/paper/Figure/SV_dataviz/ena.pdf', width=3.5, height=2.5)
drawCoolHM(ena)
dev.off()
