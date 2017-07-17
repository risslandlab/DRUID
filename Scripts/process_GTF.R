#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("One argument must be supplied (input gtf).\n", call. = FALSE)
}
if ( !grep("\\.gtf$", args[1], perl=TRUE)) {
  stop("Please ensure supplied file has suffix .gtf.\n", call. = FALSE)
}
gtf <- args[1]

require('GenomicFeatures')
require('rtracklayer')
require('plyr')

txdb <- makeTxDbFromGFF(gtf)

# Pick a transcript with the longest processed length as the representative 
# transcript for each gene.
transcripts.length <- transcriptLengths(txdb)
transcripts.length <- ddply(
  transcripts.length,
  .(gene_id),
  function(x) {
    x[which.max(x$tx_len), ]
  }
)
row.names(transcripts.length) <- transcripts.length$tx_name

# Prepare required fields for GTF file and save to disk.
exons <- exons(
  txdb,
  columns = c('gene_id', 'tx_name', 'exon_rank', 'exon_name'),
  vals = list(tx_id=transcripts.length$tx_id)
)
mcols(exons)$type <- 'exon'
mcols(exons)$source <- 'refGene'
mcols(exons)$gene_name <- mcols(exons)$gene_id
names(mcols(exons))[2] <- 'transcript_id'
names(mcols(exons))[3] <- 'exon_number'
names(mcols(exons))[4] <- 'exon_id'

export(
  exons,
  sub("\\.gtf$", ".exons.gtf", gtf, perl = TRUE),
  format='gtf'
)

# Also save a version where no exons overlap.
overlaps <- findOverlaps(
  exons,
  type = "any",
  ignoreSelf = TRUE,
  select = "arbitrary"
)
exons <- exons[is.na(ov)]

export(
  exons,
  sub("\\.gtf$", ".exons.no.gtf", gtf, perl = TRUE),
  format = 'gtf'
)

# Pick introns corresponding to representative transcript above.
introns <- intronsByTranscript(txdb, use.names = TRUE)
introns <- unlist(introns)
mcols(introns)$transcript_id <- names(introns)
mcols(introns)$exon_number <- unlist(
  sapply(
   rle(mcols(introns)$transcript_id)$lengths,
   function(x) {
     seq(1, x)
   }
  )
)
introns <- introns[mcols(introns)$transcript_id %in% transcripts.length$tx_name]

# Remove introns overlapping with any annotated exon.
exons.all <- exonsBy(txdb, by = "tx", use.names = TRUE)
overlaps <- findOverlaps(
  subject = exons.all,
  query = introns,
  select = 'first'
)
introns <- introns[is.na(ov)]

# Prepare required fields for GTF file and save to disk.
mcols(introns)$type <- 'exon'
mcols(introns)$source <- 'refGene'
mcols(introns)$gene_id <- transcripts.length[
  mcols(introns)$transcript_id, "gene_id"
]
# Rearrange columns to be consistent when writing gtf.
mcols(introns) <- mcols(introns)[, c(5, 1:4)]
mcols(introns)$exon_id <- paste(
  substr(  # transcript without suffix
    x = mcols(introns)$transcript_id,
    start = 0,
    stop = nchar(mcols(introns)$transcript_id ) - 1
  ),
  paste(
    mcols(introns)$exon_number,
    substr(  # organism suffix
      x = mcols(introns)$transcript_id,
      start = nchar(mcols(introns)$transcript_id),
      stop = nchar(mcols(introns)$transcript_id)
    ),
    sep = ''
  ),
  sep = '.'
)
mcols(introns)$gene_name <- mcols(introns)$gene_id
names(introns) <- NULL

export(
  introns,
  sub("\\.gtf$", ".introns.gtf", gtf, perl = TRUE),
  format = 'gtf'
)
