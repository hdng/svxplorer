# Summarize sample count per sliding window from bed overlap
# Created by: Ha X. Dang <haxdang attt gmail dottt com>

args = commandArgs(T)
#args = c('output/genome.win-100kb.step-1kb.bed.overlap.tsv', 'output', 'T')
if (length(args) != 3){
    stop('Rscript thisfile <overlap-file> <output-dir> <plot:T/F>\n')
}
overlap.file = args[1]
out.dir = args[2]
plotit = as.logical(args[3])
if (is.na(plotit)){plotit = F}

library(grid)
library(gridBase)
library(gridExtra)
library(ggplot2)
library(reshape2)


# read window vs break point overlap from bedtools
cat('Reading window data...\n')
x = read.table(overlap.file, header=F, sep='\t', stringsAsFactors=F)
x = x[, c(1,2,3,7)]
colnames(x) = c('chr', 'start', 'stop', 'sample')
x$svtype = gsub('^.*/', '', x$sample)
x$sample = gsub('/.*', '', x$sample)
#svtypes = c('BND', 'DEL', 'DUP', 'NS', 'INV')
svtypes = setdiff(unique(x$svtype), '.')
#x = x[x$svtype %in% c(svtypes, '.'),]

samples = setdiff(unique(x$sample), '.')
total.samples = length(samples)

#chrs = c('chr21')
#x = x[x$chr %in% chrs,]

# collect all windows (should have everything including windows w/o break points)
w = x[, c('chr', 'start', 'stop')]
w = w[!duplicated(w),]
rownames(w) = paste(w$chr, w$start, w$stop, sep='_')
cat('Found total', nrow(w), 'unique windows\n')

# count number of samples having at least one break point in each window
# return data.frame of sample count and details for all windows
countSamples <- function(x, svtype='ALL'){
    cat('Counting samples for ', svtype, '...')
    if (svtype != 'ALL'){
        z = x[x$svtype == svtype,]
    }else{
        z = x
    }
    z = z[, 1:4]
    
    # make sure each sample is reported only once for each window
    z = z[!duplicated(z),]

    z1 = aggregate(sample ~ chr+start+stop, data=z, FUN=paste, collapse=',')
    z2 = aggregate(sample ~ chr+start+stop, data=z, FUN=length)
    colnames(z2)[4] = 'num.samples'
    zz = merge(z2, z1)
    zz$num.samples[zz$sample == '.'] = 0
    if (svtype != 'ALL'){
        missing.windows = w[!(rownames(w) %in% paste(zz$chr, zz$start, zz$stop, sep='_')),]
        missing.windows$num.samples = as.integer(0)
        missing.windows$sample = '.'
        zz = rbind(zz, missing.windows)
    }

    zz$pct.samples = zz$num.samples/total.samples*100
    zz$pos = (zz$start + zz$stop)/2
    zz$svtype = svtype
    cat(nrow(zz), 'windows reported (including windows w/o count).\n')
    if(nrow(zz) != nrow(w)){stop('ERR: something wrong, number of reported windows != total windows.\n')}
    return(zz)
}

# summarize by types separately, then combine. This is slow!
cat('Summarizing...\n')
aaa = countSamples(x)
for (svtype in svtypes){
    aaa = rbind(aaa, countSamples(x, svtype))
}
zz = aaa
saveRDS(zz, file=file.path(out.dir, 'count.rds'))

#stop()

# read genes of interest to plot along with results
g = read.table('output/genes-of-interest.bed', header=F, stringsAsFactors=F)
colnames(g) = c('chr', 'start', 'stop', 'gene')
g$pos = (g$start + g$stop)/2


# read breakpoints (to plot)
cat('Readding original breakpoints...\n')
bp = read.table('output/bp.bed', header=F, stringsAsFactors=F, sep='\t')
colnames(bp) = c('chr', 'start', 'stop', 'sample', 'score')
bp$type = gsub('^.*/', '', bp$sample)
bp$sample = gsub('/.*', '', bp$sample)
bp$pos = (bp$start + bp$stop)/2
bp = bp[bp$type %in% svtypes,]


printp <- function(p){
    plot.new()
    vps = baseViewports()
    pushViewport(vps$figure)
    vp1 = plotViewport(c(0,0,0,0))
    print(p, vp=vp1)
    popViewport()
}

# Plot bk count, sample count by chr
if (plotit){
    cat('Plotting...\n')
    dir.create(out.dir)
    plot.out.dir = file.path(out.dir , 'plots')
    dir.create(plot.out.dir)
    colors = c('red', 'blue', 'green', 'purple', 'cyan', 'orange', 'black', 'black')
    names(colors) = c('ALL', svtypes)
    for (chr in unique(zz$chr)){
        cat(chr, '\n')
        g1 = g[g$chr == chr,]
        fake = data.frame(chr=chr, pos=-1000, gene='', stringsAsFactors=F)
        if (nrow(g1) == 0){g1 = fake}
        g1$y = -5

        zz1 = zz[zz$chr == chr,]

        maxx = max(zz1$stop)
        minx = min(zz1$start) - 0.025*maxx
        maxx = maxx*1.025

        png(paste0(plot.out.dir, '/', chr, '.png'), units='in', res=300, width=11, height=9)
        par(mfrow=c(length(svtypes) + 2,1))

        # histogram plot for # of break points
        p1 = (ggplot(bp[bp$chr == chr,], aes(x=pos))
            + geom_histogram(binwidth=100000)
            + theme_bw()
            + scale_x_continuous(limits=c(minx, maxx))
            + scale_y_continuous('# breakpoints', limits=c(0, 100))
            + xlab('chr position')
        )
        printp(p1)

        # sliding window plot for # of samples for each type of break points
        for (svtype in c('ALL', svtypes)){
            zz2 = zz1[zz1$svtype == svtype & zz1$pct.samples > 0, ]
            p2 = (ggplot(zz2, aes(x=pos, y=pct.samples))
                #+ geom_bar(stat='identity', fill=colors[svtype], color=colors[svtype])
                + geom_bar(stat='identity', fill=colors[svtype], size=0)
                + theme_bw()
                + geom_text(data=g1, aes(x=pos, y=y, label=gene), angle=90, hjust=1,color='blue', size=1.5)
                + scale_x_continuous(limits=c(minx, maxx))
                + xlab(NULL)
                + scale_y_continuous('% samples', limits=c(-20, 100))
                + annotate('text', x=0, y=50, label=svtype, color=colors[svtype])
            )
            printp(p2)
        }

        dev.off()
    }
}

# prepare data for circos
cat('Preparing circos track...\n')
circos.out.dir = file.path(out.dir, '4circos')
dir.create(circos.out.dir)
for (svtype in c('ALL', svtypes)){
    cat(svtype, '\n')
    zz1 = zz[zz$svtype == svtype & zz$pct.samples > 0,]
    zz1$chr = gsub('chr', 'hs', zz1$chr)
    write.table(zz1[, c('chr', 'start', 'stop', 'pct.samples')],
        file= paste0(circos.out.dir, '/pct-samples-per-window.', svtype, '.tsv'),
        row.names=F, sep='\t', quote=F)
}


