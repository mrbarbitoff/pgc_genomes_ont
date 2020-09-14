library(ggplot2)
library(ape)
library(cowplot)


ty = read.table('cumulative_Ty.tsv', sep='\t', header=F)
head(ty)
ty$V4[ty$V4 == 'Ty1,Ty2'] = 'Ty1'
ty$V2 = factor(ty$V2, levels=c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
                               'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'))

dummy = read.table('dummy.tsv', sep='\t', header=F)

ggplot(ty, aes(x=V3)) + geom_vline(aes(xintercept=V3, col=V4), lwd=1) +
  facet_grid(vars(V1), vars(V2), scales='free') + theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = 'top') + xlab('') +
  geom_blank(data=dummy, aes(x=V3))

# Trying different widths for chromosomes

p <- ggplot(ty, aes(x=V3)) + geom_vline(aes(xintercept=V3, col=V4), lwd=1) +
  facet_grid(vars(V1), vars(V2), scales='free') + theme_bw() +
  theme(axis.text.x = element_blank(), 
        legend.position = 'top',
        panel.grid = element_blank()) + xlab('') +
  geom_blank(data=dummy, aes(x=V3))

gp <- ggplotGrob(p)
facet.columns <- seq(5, 35, by = 2)

true_widths = dummy[dummy$V3 > 1, ][1:16, 'V3']
true_widths = true_widths/max(true_widths)
scaler = min(true_widths)
true_widths = true_widths/scaler

gp$widths[facet.columns] <- gp$widths[facet.columns] * true_widths
grid::grid.draw(gp)



# Trees
shared = read.table('shared_Ty.tsv', sep='\t', header=F)
head(shared)

shared_matrix = matrix(1 - shared$V3/shared$V4, nrow=4, ncol=4)
rownames(shared_matrix) = c('1A', '74D', 'S288C', 'W303')
colnames(shared_matrix) = c('1A', '74D', 'S288C', 'W303')
print(shared_matrix)

tr <- nj(as.dist(shared_matrix))
plot(tr, "u")


snps = read.table('snp_stats.tsv', sep='\t', header=F)
head(snps)

snp_matrix = matrix(snps$V3, nrow=4, ncol=4)
rownames(snp_matrix) = c('1A', '74D', 'S288C', 'W303')
colnames(snp_matrix) = c('1A', '74D', 'S288C', 'W303')
print(snp_matrix)

tr <- nj(as.dist(snp_matrix))
plot(tr, "u")

shared$snps = snps$V3/1000
shared$ty_prop = shared$V3/shared$V4
to_corr = shared[c(2, 3, 4, 7, 8, 12), ]
ggplot(to_corr, aes(x=ty_prop, y=snps)) + geom_point(size=3) +
  geom_smooth(method='lm', col='black', lwd=0.6) +
  theme_bw() + xlab('Shared Ty element locations (%)') +
  ylab('Number of substitutions (x1000)')


