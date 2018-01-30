x = read.table('GRCh38Decoy_refseq_genelocs_from_refFlat.bed', header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(x) = c('chr', 'start', 'stop', 'gene', 'score', 'strand')

x = x[x$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]

x$gene = gsub('~.*$', '', x$gene)

x1 = x[order(x$start, x$gene), c('chr', 'start', 'gene', 'strand')]
x1 = x1[!duplicated(x1$gene),]
x1 = x1[order(x1$gene),]

x2 = x[order(-x$stop, x$gene), c('chr', 'stop', 'gene', 'strand')]
x2 = x2[!duplicated(x2$gene),]
x2 = x2[order(x2$gene),]


#sanity check
table(x1$gene == x2$gene)

z = merge(x1,x2)
print(table(duplicated(z$gene)))
z$score = 1000
#refFlat by UCSC is 0-based start (similar to bed)
z = z[, c('chr', 'start', 'stop', 'gene', 'score', 'strand')]

outGeneBedFile = 'genes.bed'
write.table(z, file=outGeneBedFile, row.names=F, col.names=F, sep='\t', quote=F)

# write genes with multiple strands and chrs
diff.chr.strand.genes = x1$gene[x1$chr != x2$chr | x1$strand != x2$strand]
write.table(x[x$gene %in% diff.chr.strand.genes,],
    paste0(outGeneBedFile, '.multi-chr-or-strand.tsv'), sep='\t', quote=F, row.names=F)

