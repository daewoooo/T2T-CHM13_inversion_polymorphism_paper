## Load required libraries
library(GenomicRanges)
library(dplyr)
library(parallelDist)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(aplot)

## Load required code to load and process PAF alignments
sapply(list.files(pattern="[.]R$", path="/Porubsky_etal_2023/Rfunctions/", full.names=TRUE), source)

## Get ROI on CHM13v1.1
plt.region <- as('chr15:81700000-83500000', 'GRanges')

## Plot all alignments ##
#########################
paf.files <- list.files('/Porubsky_etal_2023/Data/15q25.2_align2T2T-CHM13v1.1/', pattern = '\\.paf$', full.names = TRUE)

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
library(parallelDist)
library(ggtree)
library(ape)
library(phangorn)

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

plt <- ggtree(tree) + geom_tiplab()
plt1 <- plt +  
  scale_x_continuous(limits = c(0, max(tree$edge.length) + 0.1)) +
  ggtitle('15q25.2 distance tree') +
  theme(legend.position = 'left')

d <- fortify(tree)
dd <- subset(d, isTip)
ordered.labels <- dd$label[order(dd$y, decreasing=TRUE)]
ordered.labels <- rev(ordered.labels)

matches.df$ID <- factor(matches.df$ID, levels = ordered.labels)
colors <- list('+' = 'chartreuse4', '-' = 'darkgoldenrod2')

plt2 <- ggplot() +
  geom_linerange(data=matches.df, aes(x=ID, ymin=start, ymax=end, color=strand), alpha=0.5, size=2) +
  scale_color_manual(values = colors) +
  scale_y_continuous(labels=comma, limits = c(start(plt.region), end(plt.region)), expand = c(0,0), name='') +
  ylab('') +
  coord_flip() +
  theme_minimal() +
  theme(#axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

final.plt <- plot_grid(plt1, plt2, nrow = 1, align = 'h', axis = 'tb')

## Make minimap dotplot [Region spanning pairs only] ##
#######################################################
inputfolder <- '/Porubsky_etal_2023/Data/15q25.2_selfAlignments/'
mm2.coords.files <- list.files(path = inputfolder, pattern = '\\.paf$', full.names = TRUE)

min.align.len <- 5000
min.selfaln.dist <- 500000 ## Distance set 500kb to ensure only SD pairs spanning ROI are considered

summary.l <- list()
for (i in seq_along(mm2.coords.files)) {
  mm2.coords <- mm2.coords.files[i]
  asm.id <- basename(mm2.coords)
  asm.id <- gsub(asm.id, pattern = '.paf', replacement = '')
  message('Processing: ', asm.id)
  
  ## Load self-alignments
  paf.table <- crunchPaf(paf.file = mm2.coords, min.align.len = min.align.len, min.selfaln.dist = min.selfaln.dist, bin.paf.aln = FALSE)

  ## Store alignments
  aln.data <- paf.table$M
  if (nrow(aln.data) > 0) {
    dir <- ifelse(aln.data$strand == '+', 'forw', 'rev')    
    summary.df <- data.frame(aln.len=aln.data$aln.len, dir=dir)
    summary.df <- summary.df %>% group_by(dir) %>% summarise(total.bp = sum(aln.len)) %>% mutate(aln.len.fraction = total.bp / sum(total.bp))
    summary.df$asm.id <- asm.id
    summary.l[[length(summary.l) + 1]] <- summary.df
  }  
}  
summary.selfaln.df <- do.call(rbind, summary.l)

## Plot self-alignment summary ##
#################################
asm.ids <- unique(summary.selfaln.df$asm.id)
sample.ord <- summary.selfaln.df %>% filter(dir == 'forw') %>% arrange(desc(aln.len.fraction))
sample.ord <- sample.ord$asm.id
sample.ord <- c(sample.ord, asm.ids[!asm.ids %in% sample.ord])
summary.selfaln.df$asm.id <- factor(summary.selfaln.df$asm.id, levels = sample.ord)
summary.selfaln.df$ID <- 'Human'

p1 <- ggplot(summary.selfaln.df) +
  geom_col(aes(y=asm.id, x=aln.len.fraction, fill=dir)) +
  #geom_point(aes(y=0, x=asm.id, color=inv.status)) +
  geom_vline(xintercept = c(0.25, 0.5, 0.75), color='white', linetype='dotted', size=1) +
  scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2'), name='Alignment\ndirection') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
  scale_x_continuous(expand = c(0, 0.01), name = 'Fraction of basepairs') +
  ylab('Assembly ID') +
  theme_minimal()

p2 <- ggplot(summary.selfaln.df) +
  geom_col(aes(y=asm.id, x=total.bp, fill=dir)) +
  #geom_point(aes(y=0, x=asm.id, color=inv.status)) +
  geom_vline(xintercept = c(0.25, 0.5, 0.75), color='white', linetype='dotted', size=1) +
  scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2'), name='Alignment\ndirection') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
  scale_x_continuous(expand = c(0, 0.01), name = 'Total # of basepairs (Stacked)') +
  ylab('Assembly ID') +
  theme_minimal()

## Merge all plots together
final.plt <- plt2 %>% insert_left(plt1, width = 1) %>% insert_right(p2, width = 0.5) %>% insert_right(p1, width = 0.5)
