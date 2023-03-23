## Load required libraries and functions
library(GenomicRanges)
library(regioneR)
library(dplyr)
library(ggplot2)
  
source("/Porubsky_etal_2023/Rfunctions/collapseBins.R")

chromosomes <- paste0('chr', c(1:22, 'X'))

## Load CHM13 centromere annotation track
cent.df <- read.table("/Porubsky_etal_2023/Data/chm13.draft_v1.1.cenAnnotation.bed", stringsAsFactors = FALSE)
cent.df <- cent.df[,c(1,2,3,4)]
colnames(cent.df) <- c('seqnames', 'start', 'end', 'censat')
## Assign CenSat category
cent.df$censat_categ <- ''
cent.df$censat_categ[grep(cent.df$censat, pattern = 'hsat')] <- 'human-satellites' 
cent.df$censat_categ[grep(cent.df$censat, pattern = 'bsat')] <- 'beta-satellites'
cent.df$censat_categ[grep(cent.df$censat, pattern = 'hor')] <- 'alpha-satellites HOR array'
cent.df$censat_categ[grep(cent.df$censat, pattern = 'mon')] <- 'alpha-satellite monomeric regions'
cent.df$censat_categ[grep(cent.df$censat, pattern = 'ct')] <- 'centric transition regions'
cent.df$censat_categ[grep(cent.df$censat, pattern = 'rDNA')] <- 'rDNA'
cent.df <- cent.df[cent.df$censat_categ != '',]
cent.df$seqnames <- factor(cent.df$seqnames, levels = paste0('chr', c(1:22, 'X')))
## Remove centric transition regions
cent.df <- cent.df[cent.df$censat_categ != 'centric transition regions',]
cent.gr <- makeGRangesFromDataFrame(cent.df, keep.extra.columns = TRUE)
names(cent.gr) <- NULL
## Get centromeric regions per chromosome
cent.region.gr <- cent.gr[cent.gr$censat_categ %in% c('human-satellites', 'alpha-satellites HOR array', 'rDNA')]
#cent.region.gr <- cent.gr[grep(cent.gr$censat, pattern = '^hor_')]
cent.region.gr$chrom <- as.character(seqnames(cent.region.gr))
cent.region.gr <- collapseBins(cent.region.gr, id.field = 3)
cent.region.gr <- cent.region.gr[,0]
cent.region.gr <- sort(cent.region.gr)
## Add 1Mbp at each side of each centromeric region
cent.region.gr <- resize(cent.region.gr, width = width(cent.region.gr) + 2000000, fix = 'center')
start(cent.region.gr) <- pmax(1, start(cent.region.gr))

## Load CHM13-T2T Strand-seq inversion callset ##
chm13.data.df <- read.table("/Porubsky_etal_2023/Data/inversionCallset_T2TCHM13v1.1coords.tsv", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
chm13.inv.gr <- makeGRangesFromDataFrame(chm13.data.df, keep.extra.columns = TRUE)
chm13.inv.gr$id <- as.character(chm13.inv.gr)
## Select balanced inversions only
chm13.inv.gr <- chm13.inv.gr[chm13.inv.gr$categ %in% c('inv')]
chm13.inv.gr <- sort(chm13.inv.gr)

## Get enrichment of inversions counts within pericentromeres using regioneR [CHM13] ##
#######################################################################################
## Get genome lengths from fasta index
fai <- '/Porubsky_etal_2023/Data/chm13_v1.1_plus38Y.fasta.fai'
fai.tab <- utils::read.table(fai)
fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
whole.genome <- data.frame(seqnames=fai.tab$V1, start=1, end=fai.tab$V2)
whole.genome.gr <- makeGRangesFromDataFrame(whole.genome)

## Calculate enrichment per chromosome
permuted.l <- list()
observed.l <- list()
pval.l <- list()
for (i in seq_along(cent.region.gr)) {
  cent.gr <- cent.region.gr[i]
  inv.gr <- keepSeqlevels(chm13.inv.gr, value = seqnames(cent.gr), pruning.mode = 'coarse')
  chromosome <- as.character(seqnames(cent.gr))
  message('Processing chromosome: ', chromosome)
  pt <- permTest(
    A=inv.gr, 
    B=cent.gr, 
    randomize.function=randomizeRegions,
    evaluate.function=numOverlaps,
    genome=whole.genome[whole.genome$seqnames == chromosome,], 
    ntimes=1000, 
    allow.overlaps=FALSE, 
    count.once=FALSE,
    mc.set.seed=FALSE,
    mc.cores=16)
  ## Store results
  perm <- data.frame(value=pt$numOverlaps$permuted)
  obs <- data.frame(value=pt$numOverlaps$observed)
  pval <- data.frame(value=pt$numOverlaps$pval)
  perm$chrom <- chromosome
  obs$chrom <- chromosome 
  pval$chrom <- chromosome 
  permuted.l[[i]] <- perm
  observed.l[[i]] <- obs
  pval.l[[i]] <- pval
}  
permuted.df <- do.call(rbind, permuted.l)
observed.df <- do.call(rbind, observed.l)
pval.df <- do.call(rbind, pval.l)
## Adjust p-value for multiple testing
pval.df$pVal_BonfAdjust <- p.adjust(pval.df$value, method = "bonferroni", n = length(pval.df$value))
pval.df$label <- round(pval.df$pVal_BonfAdjust, digits = 3)

permuted.df$chrom <- factor(permuted.df$chrom, paste0('chr', c(1:22, 'X')))
observed.df$chrom <- factor(observed.df$chrom, paste0('chr', c(1:22, 'X')))

## Get fold-enrichment for chromosomes with signif p-values
mean.vals <- permuted.df %>% group_by(chrom) %>% summarise(mean.perm = mean(value))
mean.vals$observed <- observed.df$value[match(observed.df$chrom, mean.vals$chrom)]
mean.vals$fold.change <- round(mean.vals$observed / mean.vals$mean.perm,digits = 1)
mean.vals$label <- ''
mean.vals$label[!is.nan(mean.vals$fold.change)] <- paste0(mean.vals$fold.change[!is.nan(mean.vals$fold.change)], 'x')
mean.vals$pval_label <- pval.df$label[match(pval.df$chrom, mean.vals$chrom)]

## Visualize enrichment analysis
my_theme <- theme(legend.position="right",
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  panel.background = element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  strip.text.y.left = element_text(angle = 0),
                  strip.background =element_blank(),
                  panel.spacing.y=unit(1, "mm"))

x.breaks <- 0:(max(observed.df$value) + 1)
plt <- ggplot(data=permuted.df, aes(x=value, y=1)) +
  geom_violin(fill='black') +
  geom_point(data=observed.df, color='red') +
  geom_text(data=mean.vals[mean.vals$pval_label < 0.05,], aes(x=observed, y=1, label=label), inherit.aes = FALSE, hjust=-0.5) +
  scale_x_continuous(expand = c(0,0), breaks = x.breaks, limits = c(0, max(x.breaks))) +
  ylab('') +
  xlab('Pericentromeric inversion count') +
  facet_grid(chrom ~ ., switch='y') +
  my_theme
