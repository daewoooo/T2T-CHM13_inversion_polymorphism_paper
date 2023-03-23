## Extract repeats from inverted and direct haplotypes ##
#########################################################
library(Biostrings)
library(SVbyEye)
library(DECIPHER)
library(GenomicRanges)
library(ggmsa)

## Load required code to load and process PAF alignments
sapply(list.files(pattern="[.]R$", path="/Porubsky_etal_2023/Rfunctions/", full.names=TRUE), source)

## Get inverted haplotype HG02257_1
mm2.self.aln <- '/Porubsky_etal_2023/Data/15q25.2_selfAlignments/HG02257_1.paf'

min.align.len <- 10000
min.selfaln.dist <- 500000 
inv.hap <- crunchPaf(paf.file = mm2.self.aln, min.align.len = min.align.len, min.selfaln.dist = min.selfaln.dist, bin.paf.aln = FALSE)
inv.hap <- inv.hap$M
## Keep only reverse oriented pairs
inv.hap <- inv.hap[inv.hap$strand == '-',]
gr1 <- makeGRangesFromDataFrame(inv.hap, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
gr2 <- makeGRangesFromDataFrame(inv.hap, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
gr1 <- range(gr1)
gr2 <- range(gr2)

asm.fasta <- '/Porubsky_etal_2023/Data/15q25.2_assemblies/chr15_81700000_83500000_HG02257_1.fasta'
fa.file <- open(Rsamtools::FaFile(asm.fasta))
fa.idx <- Rsamtools::scanFaIndex(fa.file)
seqlevels(gr1) <- seqlevels(fa.idx)
inv.seq1 <- Rsamtools::scanFa(file = fa.file, param = gr1, as = "DNAStringSet")
seqlevels(gr2) <- seqlevels(fa.idx)
inv.seq2 <- Rsamtools::scanFa(file = fa.file, param = gr2, as = "DNAStringSet")

## Get direct haplotype CHM13v2_1
mm2.self.aln <- '/Porubsky_etal_2023/Data/15q25.2_selfAlignments/CHM13v2_1.paf'

min.align.len <- 10000
min.selfaln.dist <- 500000 
inv.hap <- crunchPaf(paf.file = mm2.self.aln, min.align.len = min.align.len, min.selfaln.dist = min.selfaln.dist, bin.paf.aln = FALSE)
inv.hap <- inv.hap$M
## Keep only reverse oriented pairs
inv.hap <- inv.hap[inv.hap$strand == '-',]
gr1 <- makeGRangesFromDataFrame(inv.hap, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
gr2 <- makeGRangesFromDataFrame(inv.hap, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
gr1 <- range(gr1)
gr2 <- range(gr2)

asm.fasta <- '/Porubsky_etal_2023/Data/15q25.2_assemblies/chr15_81700000_83500000_CHM13v2_1.fasta'
fa.file <- open(Rsamtools::FaFile(asm.fasta))
fa.idx <- Rsamtools::scanFaIndex(fa.file)
seqlevels(gr1) <- seqlevels(fa.idx)
dir.seq1 <- Rsamtools::scanFa(file = fa.file, param = gr1, as = "DNAStringSet")
seqlevels(gr2) <- seqlevels(fa.idx)
dir.seq2 <- Rsamtools::scanFa(file = fa.file, param = gr2, as = "DNAStringSet")

## Reverse complement opposing SD pair
dir.seq2 <- reverseComplement(dir.seq2)
inv.seq2 <- reverseComplement(inv.seq2)

seqs <- c(dir.seq1, dir.seq2, inv.seq1, inv.seq2)
names(seqs) <- paste0(names(seqs), c('dir_proximal', 'dir_distal_revcomp', 'inv_proximal', 'inv_distal_revcomp'))

## Perform MSA ##
#################
DNA <- AlignSeqs(seqs)

nucfreq.m <- t(consensusMatrix(DNA, baseOnly=TRUE))
to.select <- apply(nucfreq.m, 1, function(x) x[5] < 3 & max(x) != 4)
to.select <- which(to.select == TRUE)
## Mask positions where single sample is missing a sequence
nucfreq.dir <- t(consensusMatrix(DNA[1:2], baseOnly=TRUE))
to.remove.dir <- apply(nucfreq.dir, 1, function(x) x[5] > 0)
to.remove.dir <- which(to.remove.dir == TRUE)
nucfreq.inv <- t(consensusMatrix(DNA[3:4], baseOnly=TRUE))
to.remove.inv <- apply(nucfreq.inv, 1, function(x) x[5] > 0)
to.remove.inv <- which(to.remove.inv == TRUE)
## Remove missing sequence from either direct or inverted haplotypes
to.remove <- unique(c(to.remove.dir, to.remove.inv))
to.select <- to.select[!to.select %in% to.remove]

#nucfreq.m[to.select,]
s1 <- letter(DNA[[1]], to.select)
s2 <- letter(DNA[[2]], to.select)
s3 <- letter(DNA[[3]], to.select)
s4 <- letter(DNA[[4]], to.select)

DNA.sub <- DNAStringSet(c(s1, s2, s3, s4))
names(DNA.sub) <- names(DNA)

plt.df <- tidy_msa(DNA.sub)
plt.df$PSV <- 0
for (i in seq_along(unique(plt.df$position))) {
  pos <- plt.df[plt.df$position == i,]
  alleles <- unique(pos$character)
  names(alleles) <- 1:length(alleles)
  pos$PSV <- match(pos$character, alleles)
  pos$PSV[pos$character == '-'] <- 0
  plt.df[plt.df$position == i,] <- pos
}

chunks <- 20
gen.pos <- to.select
breaks <- floor(length(gen.pos)/chunks)
breaks <- c(0, breaks*1:chunks)
labels <- c(0, gen.pos[breaks])

plt.df$PSV <- as.character(plt.df$PSV)
plt <- ggplot(plt.df) + 
  geom_tile(aes(x=position, y=name, fill=PSV)) +
  scale_x_continuous(breaks = breaks, labels = labels, expand = c(0,0), name = 'SD position', limits = c(1, max(plt.df$position))) +
  ylab('SD ids') +
  scale_fill_manual(values = c('0' = 'white', '1' = 'lightseagreen', '2' = 'yellow3', '3' = 'black', '4' = 'red'), name="Alleles (PSVs)") +
  theme_bw()
