## Load required libraries
library(GenomicRanges)
library(dplyr)
library(parallelDist)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)

## Load required code to load PAF alignments
source('/Porubsky_etal_2023/Rfunctions/loadMinimapAlignments.R')

## Get ROI on CHM13v1.1
plt.region <- as('chrX:147500000-148500000', 'GRanges')

## Load SD annotation ##
SDs.df <- read.table("/Porubsky_etal_2023/Data/CHM13-T2T_sedef_v1.1.gz")
SDs.df <- SDs.df[,c(1,2,3,24)]
colnames(SDs.df) <- c('seqnames', 'start', 'end', 'fracMatch')
SDs.gr <- makeGRangesFromDataFrame(SDs.df, keep.extra.columns = TRUE)
## Subset SD annotation ##
sd.annot.gr <- subsetByOverlaps(SDs.gr, plt.region)
sd.annot.df <- as.data.frame(sd.annot.gr)
## Define SD colors
sd.categ <- findInterval(sd.annot.df$fracMatch, vec = c(0.95, 0.98, 0.99))
sd.categ <- dplyr::recode(sd.categ, '0' = '<95%', '1' = '95-98%', '2' = '98-99%', '3'='>=99%')
sd.categ <- factor(sd.categ, levels=c('<95%', '95-98%', '98-99%', '>=99%'))
sd.annot.df$sd.categ <- sd.categ

## Plot all alignments ##
#########################
paf.files <- list.files('/Porubsky_etal_2023/Data/Xq28_align2T2T-CHM13v1.1/', pattern = '\\.paf$', full.names = TRUE)
## Use Clint PRT as outgroup
clint.paf.files <- paf.files[grep(paf.files, pattern = 'Clint')]

data.grl <- list()
for (i in seq_along(paf.files)) {
  paf.file <- paf.files[i]
  paf.id <- gsub(basename(paf.file), pattern = '.paf', replacement = '')
  message('Processing assembly: ', paf.id)
  data.gr <- cigar2ranges(paf.file = paf.file, min.insertion.size = 1000, min.deletion.size = 1000, coordinate.space = 'reference')
  ## Keep only ranges overlaping roi
  roi.gr <- GRanges(seqnames = seqnames(data.gr)[1], ranges=ranges(plt.region))
  data.gr <- subsetByOverlaps(data.gr, roi.gr)
  ## Set boundaries for plotting based on the roi.gr (plt.region)
  start(data.gr)[start(data.gr) < start(plt.region)] <- start(roi.gr) 
  end(data.gr)[end(data.gr) > end(plt.region)] <- end(roi.gr) 
  ## Resize insertion positions in reference based on their size in query
  data.gr[data.gr$cg == 'I'] <- resize(data.gr[data.gr$cg == 'I'], width = data.gr$size[data.gr$cg == 'I'], fix = 'center')
  data.gr$ID <- paf.id
  ## Get inversion status
  data.grl[[i]] <- data.gr
}
plt.gr <- do.call(c, data.grl)
matches.df <- as.data.frame(plt.gr[plt.gr$cg == 'M'])
mismatches.df <- as.data.frame(plt.gr[plt.gr$cg == 'X'])
insertions.df <- as.data.frame(plt.gr[plt.gr$cg == 'I'])
insertions.df$midpoint <- insertions.df$start + ((insertions.df$end - insertions.df$start) / 2)
aligns.gr <- range(plt.gr[plt.gr$cg == 'M'], ignore.strand=TRUE)
aligns.df <- as.data.frame(aligns.gr)

colors <- list('+' = 'chartreuse4', '-' = 'darkgoldenrod2')
plt <- ggplot() +
  geom_segment(data=aligns.df, aes(x=start, xend=end, y=0.5, yend=0.5), color='gray35') +
  geom_rect(data=matches.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=strand), alpha=0.5) +
  geom_rect(data=insertions.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+0.5), color='blue', fill='blue') +
  geom_linerange(data=insertions.df, aes(x=midpoint, ymin=level-1, ymax=level), color='blue', linetype='dotted') +
  scale_fill_manual(values = colors) +
  scale_x_continuous(labels=comma, limits = c(start(plt.region), end(plt.region)), expand = c(0,0), name='') +
  ylab('') +
  facet_grid(ID ~ ., space='free', scale='free', switch = 'y') +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Cluster haplotypes based on structural similarity ##
#######################################################
matches.gr <- plt.gr[plt.gr$cg == 'M']
matches.grl <- split(matches.gr, matches.gr$ID)

region.size <- width(plt.region)

aln.matrix <- matrix(data = 0, nrow = length(matches.grl), ncol = region.size)
for (i in seq_along(matches.grl)) {
  aln <- matches.grl[[i]]
  message('Processing alignments: ', unique(aln$ID))
  pseudo.aln.id <- rep(0, region.size)
  for (j in seq_along(aln)) {
    gr <- aln[j]
    pseudo.aln <- rep(0, region.size)
    names(pseudo.aln) <- seq(from=start(plt.region), to=end(plt.region))
    if ( as.character(strand(gr)) == '+' ) {
      pseudo.aln[names(pseudo.aln) %in% seq(from=start(gr), to=end(gr))] <- 1
    } else {
      pseudo.aln[names(pseudo.aln) %in% seq(from=start(gr), to=end(gr))] <- -1
    }
    pseudo.aln.id <- pseudo.aln.id + pseudo.aln
  }
  aln.matrix[i,] <- pseudo.aln.id
}
rownames(aln.matrix) <- names(matches.grl)

dist.m <- parDist(aln.matrix, method = 'hamming', diag = FALSE)

## Get UPGMA clustering
tree <- upgma(D = dist.m)
## Convert to hclust class
hc.obj <- as.hclust.phylo(tree)
## Obtain groups allowing for max 2% difference
groups <- cutree(hc.obj, h = 0.02)
## Set chimp as outgroup
tree <- root(tree, outgroup = "Clint-PTR_1", resolve.root = TRUE)

plt <- ggtree(tree) + geom_tiplab()
plt1 <- plt +  
  scale_x_continuous(limits = c(0, max(tree$edge.length) + 0.1)) +
  ggtitle('Xq28 distance tree') +
  theme(legend.position = 'left')

d <- fortify(tree)
dd <- subset(d, isTip)
ordered.labels <- dd$label[order(dd$y, decreasing=TRUE)]
ordered.labels <- rev(ordered.labels)

matches.df$ID <- factor(matches.df$ID, levels = rev(ordered.labels))
colors <- list('+' = 'chartreuse4', '-' = 'darkgoldenrod2')
## Highlight SD blocks
sd.highlight.df <- sd.annot.df[sd.annot.df$width > 5000,]
plt2 <- ggplot() +
  geom_rect(data=matches.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=strand), alpha=0.5) +
  geom_vline(xintercept = c(sd.highlight.df$start, sd.highlight.df$end), linetype='dashed') +
  scale_fill_manual(values = colors) +
  scale_x_continuous(labels=comma, limits = c(start(plt.region), end(plt.region)), expand = c(0,0), name='') +
  ylab('') +
  facet_grid(ID ~ ., space='free', scale='free', switch = 'y') +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

final.plt <- plot_grid(plt1, plt2, nrow = 1, align = 'h', axis = 'tb')
