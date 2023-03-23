#' Read PAF from an input file
#'
#' This function takes an PAF output file from minimap2 and loads the file along
#' with user defined set of additional alignment tags (see PAF specification).
#'
#' @param paf.file A path to a PAF file containing alignments to be loaded.
#' @param include.paf.tags Set to \code{TRUE} if all additional PAF alignment tags should be included in the output.
#' @param restrict.paf.tags Define a set of PAF tag ids (e.g. NM, cg) to be reported in the output.
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols
#' @importFrom S4Vectors lapply
#' @author David Porubsky
#'
readPaf <- function(paf.file = NULL, include.paf.tags = TRUE, restrict.paf.tags = c("NM", "cg")) {
  ## Check user input ##
  if (is.null(paf.file)) {
    stop("Path to a PAF file to load is not defined!!!")
  }
  
  ## Check if file exists and is TAB delimited
  if (file.exists(paf.file)) {
    con <- file(paf.file, "r")
    first.line <- readLines(con, n = 1)
    n.fields <- length(stringr::str_split(first.line, "\t")[[1]])
    if (n.fields < 12) {
      stop("User defined 'paf.file' has less then 12 expected tab-delimeted fields!!!")
    }
    close(con)
  } else {
    stop("User defined 'paf.file' doesn't exist!!!")
  }
  
  ## Load PAF file ##
  if (file.exists(paf.file)) {
    ## Read PAF lines
    paf.lines <- readLines(paf.file)
    if (include.paf.tags) {
      fields <- stringr::str_split(paf.lines, "\t")
    } else {
      fields <- stringr::str_split(paf.lines, "\t", 13)
    }
    paf.fields <- S4Vectors::lapply(fields, "[", 1:12)
    field.names <- c("q.name", "q.len", "q.start", "q.end", "strand", "t.name", "t.len", "t.start", "t.end", "n.match", "aln.len", "mapq")
    
    for (i in seq_along(paf.fields)) {
      attr(paf.fields[[i]], "names") <- field.names
    }
    paf <- dplyr::bind_rows(paf.fields)
    cols.num <- c(2, 3, 4, 7:12)
    paf[cols.num] <- suppressWarnings(dplyr::bind_cols(S4Vectors::lapply(paf[cols.num], as.numeric)))
    
    if (include.paf.tags) {
      if (any(lengths(fields) > 12)) {
        paf.tags <- S4Vectors::lapply(fields, function(x) paste(x[13:length(x)]))
        paf <- dplyr::bind_cols(paf, processPafTags(paf.tags = paf.tags, restrict.paf.tags = restrict.paf.tags))
      }
    }
    return(paf)
  } else {
    stop(paste0("PAF file ", paf.file, " doesn't exists !!!"))
    return(NULL)
  }
}

#' Process PAF specific alignment tags.
#'
#' @param paf.tags A \code{list} of PAF specific tag extracted from each alignment.
#' @inheritParams readPaf
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#'
processPafTags <- function(paf.tags, restrict.paf.tags = c("NM", "cg")) {
  ## Make sure restrict.paf.tags takes only expected values
  allowed.tags <- c("tp", "cm", "s1", "s2", "NM", "MD", "AS", "SA", "ms", "nn", "ts", "cg", "cs", "dv", "de", "rl")
  restrict.paf.tags <- restrict.paf.tags[restrict.paf.tags %in% allowed.tags]
  if (length(restrict.paf.tags) == 0) {
    message(paste0("Submitted 'restrict.paf.tags' are not present in the allowed set of PAF tags: ", paste(allowed.tags, collapse = "; ")))
    restrict.paf.tags <- NULL
  }
  
  tags <- character()
  to.numeric <- integer()
  res <- list()
  n <- length(paf.tags)
  t.idx <- 0
  for (ali.idx in seq_along(paf.tags)) {
    split.tags <- stringr::str_split(paf.tags[[ali.idx]], ":")
    for (tag in split.tags) {
      if (!is.null(restrict.paf.tags) & tag[1] %in% restrict.paf.tags) {
        if (!(tag[1] %in% tags)) {
          t.idx <- t.idx + 1
          if (tag[2] %in% c("f", "H", "i")) {
            to.numeric <- c(to.numeric, t.idx)
          }
          res[[tag[1]]] <- rep(NA, n)
          tags <- c(tags, tag[1])
        }
        res[[tag[1]]][ali.idx] <- tag[3]
      }
    }
  }
  
  for (i in to.numeric) {
    res[[i]] <- as.numeric(res[[i]])
  }
  dplyr::bind_rows(res)
}

#' Function to parse CIGAR string into a set interval ranges.
#'
#' @param cigar.str A character string containing alignment represented as a CIGAR string.
#' @param coordinate.space A used defined coordinate space given CIGAR should be parsed against, either 'reference' or 'query'.
#' @importFrom GenomicAlignments explodeCigarOps explodeCigarOpLengths cigarRangesAlongReferenceSpace cigarRangesAlongQuerySpace
#' @importFrom dplyr recode
#' @return A \code{\link{IRangesList}} object.
#' @author David Porubsky
#'
parseCigarString <- function(cigar.str = NULL, coordinate.space = "reference") {
  ## Parse CIGAR string ##
  cigar.id <- GenomicAlignments::explodeCigarOps(cigar.str)[[1]]
  cigar.len <- GenomicAlignments::explodeCigarOpLengths(cigar.str)[[1]]
  ## Get Reference space coords
  if (coordinate.space == "reference") {
    cigar.ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar.str)[[1]]
    # start(cigar.ranges[cigar.id == 'I']) <- end(cigar.ranges[cigar.id == 'I'])
  } else if (coordinate.space == "query") {
    cigar.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar.str)[[1]]
  } else {
    stop("Parameter 'coordinate.space' can only take values 'reference' or 'query' !!!")
  }
  ## Translate cigar symbols
  cigar.id.trans <- dplyr::recode(cigar.id, "=" = "match", "I" = "insertion", "D" = "deletion", "X" = "mismatch")
  ## Export ranges
  irl <- split(cigar.ranges, cigar.id.trans)
  return(irl)
}


#' Function to load CIGAR string reported in PAF alignments into a set of genomic ranges.
#'
#' @param paf.file A character string containing alignment represented as a CIGAR string.
#' @param min.insertion.size A minimum size (in base pairs) of an insertion to be retained.
#' @param min.deletion.size A minimum size (in base pairs) of a deletion to be retained.
#' @param collapse.mismatches Set to \code{TRUE} if mismatches should be collapsed in order expand matched regions.
#' @inheritParams readPaf
#' @inheritParams parseCigarString
#' @importFrom GenomicRanges GRanges width shift reduce strand
#' @return A \code{\link{GRanges}} object.
#' @author David Porubsky
#'
cigar2ranges <- function(paf.file = NULL, coordinate.space = "reference", min.insertion.size = 50, min.deletion.size = 50, collapse.mismatches = TRUE) {
  ## Read in coordinates from minimap2 output in PAF format
  paf.data <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
  qname <- paste(unique(paf.data$q.name), collapse = ";")
  ## Order by target sequence
  # paf.data <- paf.data[order(paf.data$t.start),]
  
  ## Process alignments ##
  matches <- list()
  mismatches <- list()
  insertions <- list()
  deletions <- list()
  for (i in seq_len(nrow(paf.data))) {
    paf.aln <- paf.data[i, ]
    ## Parse CIGAR string ##
    cg.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = coordinate.space)
    ## Get cigar ranges and offset ranges based on alignment starting position
    if (length(cg.ranges$match) > 0) {
      match.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$match)
    } else {
      match.gr <- GenomicRanges::GRanges()
    }
    if (length(cg.ranges$mismatch) > 0) {
      mismatch.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$mismatch)
      mismatch.gr$size <- 1
    } else {
      mismatch.gr <- GenomicRanges::GRanges()
    }
    if (length(cg.ranges$deletion) > 0) {
      del.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$deletion)
      ## In case of query coordinates get deletion size from the reference coordinates
      if (coordinate.space == "query") {
        ref.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = "reference")
        del.gr$size <- GenomicRanges::width(ref.ranges$deletion)
      } else {
        del.gr$size <- GenomicRanges::width(del.gr)
      }
      ## Filter deletions by size [in reference coordinates]
      del2reduce <- del.gr[del.gr$size < min.deletion.size]
      del.gr <- del.gr[del.gr$size >= min.deletion.size]
    } else {
      del2reduce <- GenomicRanges::GRanges()
      del.gr <- GenomicRanges::GRanges()
    }
    if (length(cg.ranges$insertion) > 0) {
      ins.gr <- GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$insertion)
      ## In case of reference coordinates get insertion size from the query coordinates
      if (coordinate.space == "reference") {
        qry.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = "query")
        ins.gr$size <- GenomicRanges::width(qry.ranges$insertion)
      } else {
        ins.gr$size <- GenomicRanges::width(ins.gr)
      }
      ## Filter insertions by size [in query coordinates]
      ins2reduce <- ins.gr[ins.gr$size < min.insertion.size]
      ins.gr <- ins.gr[ins.gr$size >= min.insertion.size]
    } else {
      ins2reduce <- GenomicRanges::GRanges()
      ins.gr <- GenomicRanges::GRanges()
    }
    
    ## Collapse simple mismatches
    if (collapse.mismatches) {
      match.gr <- GenomicRanges::reduce(c(match.gr, mismatch.gr))
    } else {
      match.gr <- GenomicRanges::reduce(match.gr)
    }
    ## Collapse filtered deletions
    if (!is.null(del2reduce) & length(del2reduce) > 0) {
      match.gr <- GenomicRanges::reduce(c(match.gr, del2reduce))
    }
    ## Collapse filtered insertions
    if (!is.null(ins2reduce) & length(ins2reduce) > 0) {
      match.gr <- GenomicRanges::reduce(c(match.gr, ins2reduce))
    }
    ## Set size of the final matched bases
    match.gr$size <- GenomicRanges::width(match.gr)
    
    ## Convert to query or target coordinates
    if (coordinate.space == "reference") {
      match.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$t.start)
      mismatch.gr <- GenomicRanges::shift(mismatch.gr, shift = paf.aln$t.start)
      del.gr <- GenomicRanges::shift(del.gr, shift = paf.aln$t.start)
      ins.gr <- GenomicRanges::shift(ins.gr, shift = paf.aln$t.start)
    } else {
      match.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$q.start)
      mismatch.gr <- GenomicRanges::shift(mismatch.gr, shift = paf.aln$q.start)
      del.gr <- GenomicRanges::shift(del.gr, shift = paf.aln$q.start)
      ins.gr <- GenomicRanges::shift(ins.gr, shift = paf.aln$q.start)
    }
    
    ## Prepare data for export
    if (length(match.gr) > 0) {
      GenomicRanges::strand(match.gr) <- paf.aln$strand
      match.gr$aln.id <- paste0("aln", i)
      match.gr$cg <- "M"
    }
    if (length(mismatch.gr) > 0) {
      GenomicRanges::strand(mismatch.gr) <- paf.aln$strand
      mismatch.gr$aln.id <- paste0("aln", i)
      mismatch.gr$cg <- "X"
    }
    if (length(ins.gr) > 0) {
      GenomicRanges::strand(ins.gr) <- paf.aln$strand
      ins.gr$aln.id <- paste0("aln", i)
      ins.gr$cg <- "I"
    }
    if (length(del.gr) > 0) {
      GenomicRanges::strand(del.gr) <- paf.aln$strand
      del.gr$aln.id <- paste0("aln", i)
      del.gr$cg <- "D"
    }
    
    matches[[length(matches) + 1]] <- match.gr
    mismatches[[length(mismatches) + 1]] <- mismatch.gr
    insertions[[length(insertions) + 1]] <- ins.gr
    deletions[[length(deletions) + 1]] <- del.gr
  }
  matches.gr <- do.call(c, matches)
  mismatches.gr <- do.call(c, mismatches)
  insertions.gr <- do.call(c, insertions)
  deletions.gr <- do.call(c, deletions)
  
  ## Add level to each separate alignment
  n.aln <- length(unique(matches.gr$aln.id))
  matches.gr$level <- rep(rep(c(1:2), times = n.aln)[seq_len(n.aln)], times = rle(matches.gr$aln.id)$lengths)
  insertions.gr$level <- matches.gr$level[match(insertions.gr$aln.id, matches.gr$aln.id)]
  deletions.gr$level <- matches.gr$level[match(deletions.gr$aln.id, matches.gr$aln.id)]
  mismatches.gr$level <- matches.gr$level[match(mismatches.gr$aln.id, matches.gr$aln.id)]
  
  ## Return
  final.gr <- c(matches.gr, mismatches.gr, insertions.gr, deletions.gr)
  return(final.gr)
}