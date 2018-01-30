###### Pipeline to discover genes recurrently affected by break points
# Created by: Ha X. Dang <haxdang att gmail dott com>
#
# This pipeline has two steps:
# 1-6. Discover and annotate regions that have aggrevated frequency of break
#      points, defined as regions having high number of samples having at least
#      one break point, compared with surrounding regions. 
#
#      Annotate regions with genes, and recount the number of samples having
#      break points for flanked regions of the genes
#
# 7- . Reinvestigate genes identified within aggrevated regions of breakpoints
#      by re-counting the number of samples within specific bases up/down stream
#      Various flanking lengths can be used
#
# INPUT: bedpe files for individual samples of the same cohort in a directory
#        named "bedpe", as follows:
#        
#        ./bedpe
#        --- sample-1.bedpe
#        --- sample-2.bedpe
#        --- sample-3.bedpe
#        ...................

### 1. Collect, filter, prepare break points from SV calls
# merge SV bedpe from all samples to a single file, take whatever Manta flags as PASS
mkdir output; for x in `ls bedpe/ | grep bedpe$ | sed 's/.bedpe//'`; do export PATIENT=$x; cat bedpe/${x}.bedpe | grep -v '#' | grep Manta | perl -ane '$pt=$ENV{PATIENT}; $id="$pt/$F[10]"; if ($F[11] eq "PASS"){print join("\t", $F[0], $F[1], $F[2], $F[3], $F[4], $F[5], $id, $F[7], $F[8], $F[9]), "\n"}'; done > output/all.bedpe

# 1a. filter SV if a breakpoint is in blacklist regions
pairToBed -type neither -a output/all.bedpe -b annot/hg38/longranger-blacklist.bed > output/all.black-list-filtered.bedpe

# 1b. filter SV if both mated breakpoints are in paired segdup regions defined in segdups.bedpe
# this is not yet comprehensive as segdups.bedpe may not have all pairs of segdups
# TODO: make sure segdups.bedpe has all pairs of segdups, or create custom script for filter
pairToPair -type notboth -is -a output/all.black-list-filtered.bedpe -b annot/hg38/longranger-segdups.bedpe > output/all.bedpe; rm -f output/all.black-list-filtered.bedpe

# 1c. collect all break points in bed format for easy overlapping
# duplication within sample may exist but should not affect
# total number of samples in later counting procedure
# TODO: remove this step and overlap directly with bedpe of SV bp pairs to count events
cat output/all.bedpe | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$7,$8,$9; print $4,$5,$6,$7,$8,$10}' > output/bp.bed

# 1d. extract break points by SV types
for type in `cut -f4 output/bp.bed | cut -f2 -d'/' | sort | uniq`; do cat output/bp.bed | grep "/$type" > output/bp.$type.bed; done

### 2. Generate sliding windows from genome
# 2a. for peak finding, small step, moderate window size
Rscript3.1.2 genome-to-sliding-window.r annot/hg38/chromsize.tsv 100000 1000 output/genome.win-100kb.step-1kb.bed

# 2b. for visualization, larger step and window, less data points
Rscript3.1.2 genome-to-sliding-window.r annot/hg38/chromsize.tsv 3000000 1000000 output/genome.win-3mb.step-1mb.bed

### 3. Overlap break points with sliding windows
intersectBed -wao -a output/genome.win-100kb.step-1kb.bed -b output/bp.bed > output/genome.win-100kb.step-1kb.bed.overlap.tsv

### 4. Summarize/count/plot number of samples having at least one break point in window
# 4a. prepare genes of interest bed file
cat annot/genes-of-interest.txt | ./select-rows-v2.pl 0 annot/hg38/genes.bed 3 > output/genes-of-interest.bed

# 4b. summarize/plot, a count.rds file will be placed in the output dir that will be used to find peak
Rscript3.1.2 summarize-sample-count.r output/genome.win-100kb.step-1kb.bed.overlap.tsv output T

### 5. Find peaks
Rscript3.1.2 find-peaks.r output/count.rds ALL ALL 100 5 output/peaks

### 6. Annotate grouped peaks by overlapping with genes flanked up/downstream by 2kb
# 6a. collect all peak groups, calling this new peaks
cat output/peaks/*.peak.group.bed | grep -v pct > output/peaks/peaks.bed
# 6b. flank genes by +/- 2kb, assuming flank regions are still within chrom
# TODO: make sure coordinates won't fall out of chr limit
cat annot/hg38/genes.bed | awk '{OFS="\t"}{$2=$2-2000;$3=$3+2000; print}' > output/genes.flank-2k.bed
# 6c. overlap flanked genes and grouped peaks
intersectBed -wo -a output/peaks/peaks.bed -b output/genes.flank-2k.bed > output/peaks-vs-genes.overlap.tsv
# 6d. make bed for genes found (if a gene overlaps multiple peaks, multiple bed lines are created)
# cat output/peaks-vs-genes.overlap.tsv | awk '{OFS="\t"}{start=$7+2000; stop=$8-2000; print $6,start,stop,$9,$5,$11}' > output/genes.overlap-peaks.bed
cat output/peaks-vs-genes.overlap.tsv | cut -f9 | ./select-rows-v2.pl 0 annot/hg38/genes.bed 3 > output/genes.overlap-peaks.bed

### 7. Re-overlap genes with break points and re-count samples using various 
# Steps 1-6 allow identification of genes near/overlapping regions of aggrevated break points
# Now we need to count EXACTLY how many samples having break points falling within some distance
# from the genes
#
# 7a. Extend the candidate genes by various distance from 1kb - 100kb (flanked canidate genes loci)
for k in 1000 2000 5000 10000 50000 100000; do export KKK=$k; cat output/genes.overlap-peaks.bed | perl -ane '$k=$ENV{KKK}; $F[1]=$F[1]-$k; $F[2]=$F[2]+$k; print join("\t", @F), "\n"' > output/genes.overlap-peaks.$k.bed; done

# 7b. Overlap flanked candidate genes' loci with break points again
for k in 1000 2000 5000 10000 50000 100000; do intersectBed -wo -a output/genes.overlap-peaks.$k.bed -b output/bp.bed | cut -f4,7-10 | sed 's/\//\t/' > output/genes.overlap-peaks.$k.overlap-bp.tsv; done

# 7c. Count num of samples having at least 1 breakpoint overlapping with flanked gene loci
for k in 1000 2000 5000 10000 50000 100000; do cat output/genes.overlap-peaks.$k.overlap-bp.tsv | cut -f1,5 | sort | uniq | cut -f1 | sort | uniq -c | sort -nr | awk '{print $2"\t"$1}' > output/genes.overlap-peaks.$k.overlap-bp.num-samples.tsv; done

### 9. Identify blacklist genes and filtering candidate genes
# 9a. blacklist genes = flagged genes + genes overlapping blacklist region within 50kb
( echo -e "TTC28\nTRNAA-AGC"; cat annot/blacklisted_genes.txt; intersectBed -wo -a output/genes.overlap-peaks.50000.bed -b annot/hg38/longranger-blacklist.bed | cut -f4 | sort | uniq | grep -Pv "^PTEN$" | grep -Pv "^ATM$" ) > output/genes.blacklist.txt
# 9b. check if any cosmic genes are blacklisted.
cat annot/cosmic-census-genes.tsv |  ./select-rows-v2.pl 0 output/genes.blacklist.txt 0 > output/genes.blacklist.in-cosmic.txt

# 9c. Filter blacklist genes from candidate list
for k in 1000 2000 5000 10000 50000 100000; do cat output/genes.blacklist.txt | ./select-rows-v2.pl 0 output/genes.overlap-peaks.$k.overlap-bp.num-samples.tsv 0 1 > output/genes.overlap-peaks.$k.overlap-bp.num-samples.filtered.tsv; cat output/genes.overlap-peaks.$k.overlap-bp.num-samples.filtered.tsv | ./select-rows-v2.pl 0 output/genes.overlap-peaks.$k.overlap-bp.tsv 0 > output/genes.overlap-peaks.$k.overlap-bp.filtered.tsv; done

# 9d. recount n of samples
for k in 1000 2000 10000 50000 100000; do ( echo -e "gene\tsample\ttype\tcount"; cat output/genes.overlap-peaks.$k.overlap-bp.filtered.tsv | cut -f1,5,6 | grep -v type | sort | uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$1}' ) > output/genes.overlap-peaks.$k.overlap-bp.filtered.count.tsv; done


### 10. Correlation analysis
# correl-analysis.r

### 11. Filter and plot candidate list
# plot.r
#
#
