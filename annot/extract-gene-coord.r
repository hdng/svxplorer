args = commandArgs(T)

#args = c('hg19/refFlat.txt', 'hg19/genes.bed')

inRefFlatFile = args[1]
outGeneBedFile = args[2]


x = read.table(inRefFlatFile, header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(x) = c('gene', 'transcript', 'chr', 'strand', 'txStart', 'txStop', 'cdsStart', 'cdsStop', 'exonCount', 'exonStarts', 'exonStops')
x$start = x$txStart
x$stop = x$txStop
x$score = 1000

x = x[x$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]

# take the smallest start pos among all transcripts from same gene
x1 = x[order(x$start, x$gene), c('gene', 'chr', 'start', 'strand')]
x1 = x1[!duplicated(x1$gene),]

# take the largest end pos among all transcripts from same gene
x2 = x[order(-x$stop, x$gene), c('gene', 'chr', 'stop', 'strand')]
x2 = x2[!duplicated(x2$gene),]

x1 = x1[order(x1$gene),]
x2 = x2[order(x2$gene),]

#sanity check
table(x1$gene == x2$gene)

z = merge(x1,x2)
print(table(duplicated(x$gene)))
z$score = 1000
#refFlat by UCSC is 0-based start (similar to bed)
z = z[, c('chr', 'start', 'stop', 'gene', 'score', 'strand')]
write.table(z, file=outGeneBedFile, row.names=F, col.names=F, sep='\t', quote=F)

# write genes with multiple strands and chrs
diff.chr.strand.genes = x1$gene[x1$chr != x2$chr | x1$strand != x2$strand]
write.table(x[x$gene %in% diff.chr.strand.genes,],
    paste0(outGeneBedFile, '.multi-chr-or-strand.tsv'), sep='\t', quote=F, row.names=F)



