## Detect likely minor alleles and misorients ##
################################################
## Load required libraries
library(GenomicRanges)
library(dplyr)

## Load data frames with Watson (wReads, minus) and Crick (cReads, plus) count per inverted region 
inv.stat.hg38.df <- get(load('/Porubsky_etal_2023/Data/hg38_directionalReads_INVandMISO.RData'))
inv.stat.chm13.df <- get(load('/Porubsky_etal_2023/Data/chm13_directionalReads_INVandMISO.RData'))

## Load GRCh38 inversion callset mapped to T2T-CHM13 reference
inv.mapped.gr <- get(load('/Porubsky_etal_2023/Data/hg38InversionCallset_mapped2CHM13.RData'))
## Remove overlapping inversions calls
## Remove two nested calls within large misorientation on chromosome 16
inv.mapped.gr <- inv.mapped.gr[!inv.mapped.gr$id %in% c('chr16:35157740-35523366', 'chr16:35226544-35592170')]
## Remove nested call on chromosome 17
inv.mapped.gr <- inv.mapped.gr[!inv.mapped.gr$id %in% c('chr17:43263148-43322883')]
## Remove largely overlapping call on chromosome 3
inv.mapped.gr <- inv.mapped.gr[!inv.mapped.gr$id %in% c('chr3:195680054-195724206')]

hg38.start <- as.numeric(sapply(inv.mapped.gr$id, function(x) strsplit(x, ':|-')[[1]][2]))
hg38.end <- as.numeric(sapply(inv.mapped.gr$id, function(x) strsplit(x, ':|-')[[1]][3]))
chm13.width <- width(inv.mapped.gr)
hg38.width <- hg38.end - hg38.start
size.diff <- abs(hg38.width - chm13.width)
## Get fraction of original size difference after mapping to CHM13
hg38.sizeFrac.diff <- round((size.diff / hg38.width) * 100, 3)

## Merge data
data.df <- rbind(inv.stat.hg38.df, inv.stat.chm13.df)
## Remove sites not mapped to CHM13
mapped.ids <- unique(data.df$region.id[data.df$region.id %in% inv.mapped.gr$id])
data.df <- data.df[data.df$region.id %in% mapped.ids,]
## Remove chromosome Y ranges
data.df <- data.df[grep(data.df$region.id, pattern = 'chrY', invert = TRUE),]
## Assign region category
data.df$categ <- inv.mapped.gr$categ[match(data.df$region.id, inv.mapped.gr$id)]

## Get GRCh38 minor alleles ##
summary.df <- data.df %>% 
  dplyr::filter(categ == 'inv') %>%
  group_by(region.id, asm.id) %>% 
  summarise(cReads.sum = sum(cReads), wReads.sum = sum(wReads), .groups = 'drop') %>%
  mutate(cReads.perc = cReads.sum / (cReads.sum + wReads.sum),
         wReads.perc = wReads.sum / (cReads.sum + wReads.sum))

summary.l <- split(summary.df, summary.df$region.id)

minor.alleles <- list()
for (i in seq_along(summary.l)) {
  inv.df <- summary.l[[i]]
  diff.crick <- ( abs(diff(inv.df$cReads.perc)) / (sum(inv.df$cReads.perc) / 2) ) * 100
  diff.watson <- ( abs(diff(inv.df$wReads.perc)) / (sum(inv.df$wReads.perc) / 2) ) * 100
          
  ## Minor allele is defined as region where there is at least 25% difference in a fraction of watson and crick reads 
  ## between hg38 and chm13 while chm13 contains higher fraction of crick reads
  if (diff.crick >= 25 & diff.watson >= 25 & inv.df$cReads.perc[inv.df$asm.id == 'chm13'] > inv.df$cReads.perc[inv.df$asm.id == 'hg38']) {
    ## Check if read ratio is similar as well
    chm13.counts <- sort(unlist(inv.df[inv.df$asm.id == 'chm13', c('cReads.sum', 'wReads.sum')]))
    hg38.counts <- sort(unlist(inv.df[inv.df$asm.id == 'hg38', c('cReads.sum', 'wReads.sum')]))
    ## Keep only regions with minimum of 20 mapped reads
    if (sum(chm13.counts) >= 20 & sum(hg38.counts) >= 20) {
      chm13.ratio <- (chm13.counts[1] + 0.1) / (chm13.counts[2] + 0.1)
      hg38.ratio <- (hg38.counts[1] + 0.1) / (hg38.counts[2] + 0.1)
      ## Get ratio similarity
      ratio.diff <- ( abs(chm13.ratio - hg38.ratio) / hg38.ratio ) * 100
      if (ratio.diff <= 25) {
        minor.alleles[[length(minor.alleles) + 1]] <- inv.df
      } 
    }  
  }       
}
minor.alleles.hg38.df <- do.call(rbind, minor.alleles)
minor.alleles.hg38.df$categ <- 'GRCh38.minor'

## Get CHM13 minor alleles ##
minor.alleles <- list()
for (i in seq_along(summary.l)) {
  inv.df <- summary.l[[i]]
  diff.crick <- ( abs(diff(inv.df$cReads.perc)) / (sum(inv.df$cReads.perc) / 2) ) * 100
  diff.watson <- ( abs(diff(inv.df$wReads.perc)) / (sum(inv.df$wReads.perc) / 2) ) * 100
  
  ## Minor allele is defined as region where there is at least 25% difference in a fraction of watson and crick reads 
  ## between hg38 and chm13 while chm13 contains higher fraction of crick reads
  if (diff.crick >= 25 & diff.watson >= 25 & inv.df$cReads.perc[inv.df$asm.id == 'hg38'] > inv.df$cReads.perc[inv.df$asm.id == 'chm13']) {
    ## Check if read ratio is similar as well
    chm13.counts <- sort(unlist(inv.df[inv.df$asm.id == 'chm13', c('cReads.sum', 'wReads.sum')]))
    hg38.counts <- sort(unlist(inv.df[inv.df$asm.id == 'hg38', c('cReads.sum', 'wReads.sum')]))
    ## Keep only regions with minimum of 20 mapped reads
    if (sum(chm13.counts) >= 20 & sum(hg38.counts) >= 20) {
      chm13.ratio <- (chm13.counts[1] + 0.1) / (chm13.counts[2] + 0.1)
      hg38.ratio <- (hg38.counts[1] + 0.1) / (hg38.counts[2] + 0.1)
      ## Get ratio similarity
      ratio.diff <- ( abs(chm13.ratio - hg38.ratio) / hg38.ratio ) * 100
      if (ratio.diff <= 25) {
        minor.alleles[[length(minor.alleles) + 1]] <- inv.df
      } 
    }  
  }       
}
minor.alleles.chm13.df <- do.call(rbind, minor.alleles)
minor.alleles.chm13.df$categ <- 'CHM13.minor'

## Get GRCh38 misorients ##
summary.df <- data.df %>%
  dplyr::filter(categ == 'miso') %>%
  group_by(region.id, asm.id) %>% 
  summarise(cReads.sum = sum(cReads), wReads.sum = sum(wReads), .groups = 'drop') %>%
  mutate(cReads.perc = cReads.sum / (cReads.sum + wReads.sum),
         wReads.perc = wReads.sum / (cReads.sum + wReads.sum))
summary.df$cReads.perc[is.nan(summary.df$cReads.perc)] <- 0
summary.df$wReads.perc[is.nan(summary.df$wReads.perc)] <- 0

summary.l <- split(summary.df, summary.df$region.id)

miso.alleles <- list()
for (i in seq_along(summary.l)) {
  inv.df <- summary.l[[i]]
  diff.crick <- ( abs(diff(inv.df$cReads.perc)) / (sum(inv.df$cReads.perc) / 2) ) * 100
  diff.watson <- ( abs(diff(inv.df$wReads.perc)) / (sum(inv.df$wReads.perc) / 2) ) * 100
  
  ## Minor allele is defined as region where there is at least 25% difference in a fraction of watson and crick reads 
  ## between hg38 and chm13 while chm13 contains higher fraction of crick reads
  if (diff.crick >= 25 & diff.watson >= 25 & inv.df$cReads.perc[inv.df$asm.id == 'chm13'] > inv.df$cReads.perc[inv.df$asm.id == 'hg38']) {
    ## Check if read ratio is similar as well
    chm13.counts <- sort(unlist(inv.df[inv.df$asm.id == 'chm13', c('cReads.sum', 'wReads.sum')]))
    hg38.counts <- sort(unlist(inv.df[inv.df$asm.id == 'hg38', c('cReads.sum', 'wReads.sum')]))
    ## Keep only regions with minimum of 20 mapped reads
    if (sum(chm13.counts) >= 20 & sum(hg38.counts) >= 20) {
      miso.alleles[[length(miso.alleles) + 1]] <- inv.df
    }   
  }       
}
miso.alleles.hg38.df <- do.call(rbind, miso.alleles)
miso.alleles.hg38.df$categ <- 'GRCh38.miso'

## Merge all results
## In this table each inverted region (region.id) has reported fraction of Watson and Crick reads in respected to 
## both GRCh38 (hg38) and T2T-CHM13v1.1 (chm13) and is classified as minor allele in hg38 or chm13, or misorient in hg38.
results.df <- rbind(minor.alleles.hg38.df, minor.alleles.chm13.df, miso.alleles.hg38.df)
head(results.df)
