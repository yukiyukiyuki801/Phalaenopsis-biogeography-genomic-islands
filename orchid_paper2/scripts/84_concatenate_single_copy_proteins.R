#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("usage: 84_concatenate_single_copy_proteins.R TRIMMED_DIR OUT_DIR")
trim_dir <- normalizePath(args[[1]])
out_dir <- args[[2]]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
taxa <- c("PAM","PAP","PEQ","DCAT","DNOB","PGUA","VPLA","APO","ASPA","RICE")

read_fasta <- function(f) {
  x <- readLines(f, warn=FALSE)
  h <- grep("^>", x)
  if (!length(h)) stop("No FASTA records in ", f)
  ends <- c(h[-1]-1L, length(x))
  nm <- sub("^>", "", x[h])
  seqs <- vapply(seq_along(h), function(i)
    paste(x[(h[i]+1L):ends[i]], collapse=""), character(1))
  tx <- sub("[| ].*$", "", nm)
  # Genome-derived proteins were prefixed with their taxon codes. The two
  # transcriptome-derived inputs retain their native TransDecoder/TSA IDs.
  tx[grepl("^NODE_", nm)] <- "PAM"
  tx[grepl("^GENW[0-9]", nm)] <- "PAP"
  names(seqs) <- tx
  seqs
}

files <- sort(list.files(trim_dir, pattern="[.]trim[.]fa$", full.names=TRUE))
if (!length(files)) stop("No trimmed alignments in ", trim_dir)
concat <- setNames(rep("", length(taxa)), taxa)
qc <- data.frame()
partitions <- data.frame()
cursor <- 1L
for (f in files) {
  s <- read_fasta(f)
  lens <- nchar(s)
  og <- sub("[.]trim[.]fa$", "", basename(f))
  reason <- "retained"
  if (!setequal(names(s), taxa)) reason <- "not_exactly_one_sequence_per_taxon"
  if (length(unique(lens)) != 1L) reason <- "unaligned_lengths"
  aln_len <- if (length(lens)) lens[[1]] else 0L
  gap_fraction <- if (aln_len > 0) mean(strsplit(paste(s,collapse=""),"")[[1]] %in% c("-","?","X")) else 1
  if (aln_len < 50L) reason <- "trimmed_length_lt50aa"
  if (gap_fraction > 0.50) reason <- "gap_or_unknown_fraction_gt0p50"
  qc <- rbind(qc, data.frame(orthogroup=og, aligned_length=aln_len,
                             gap_unknown_fraction=gap_fraction, status=reason))
  if (reason != "retained") next
  s <- s[taxa]
  concat[taxa] <- paste0(concat[taxa], s[taxa])
  partitions <- rbind(partitions, data.frame(orthogroup=og, start=cursor,
                                              end=cursor+aln_len-1L, length=aln_len))
  cursor <- cursor + aln_len
}
if (!nrow(partitions)) stop("No alignments passed QC")
write.table(qc, file.path(out_dir,"single_copy_alignment_qc.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(partitions, file.path(out_dir,"single_copy_partitions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

fa <- file.path(out_dir,"strict_single_copy_concat_aa.fasta")
con <- file(fa,"w"); on.exit(close(con),add=TRUE)
for (tx in taxa) writeLines(c(paste0(">",tx), concat[[tx]]), con)
close(con); on.exit(NULL,add=FALSE)

phy <- file.path(out_dir,"strict_single_copy_concat_aa.phy")
con <- file(phy,"w"); on.exit(close(con),add=TRUE)
writeLines(sprintf("%d %d",length(taxa),nchar(concat[[1]])),con)
for (tx in taxa) writeLines(sprintf("%-10s %s",tx,concat[[tx]]),con)
close(con); on.exit(NULL,add=FALSE)

write.table(data.frame(n_input_alignments=length(files), n_retained=nrow(partitions),
                       n_taxa=length(taxa), concatenated_aa=nchar(concat[[1]])),
            file.path(out_dir,"strict_single_copy_concat_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
