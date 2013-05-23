##----------------------------------------------------------------------------
##  Data processing functions for VDJ sequencing data from Adaptive
##  Erik Clarke (ecl@mail.med.upenn.edu)
##----------------------------------------------------------------------------
library(RMySQL)
library(seqinr)
library(Biostrings)

## Constants
##----------------------------------------------------------------------------

# Load MySQL connection info
# This should have variables named my.user, my.pass, my.host, my.db, my.table
source("mysql_info.R")

# File info
root <- "~/vdj_xscid"
unfiltered.path <- "data/patients/processed/unfiltered"
filtered.path <- "data/patients/processed/filtered"
controls.path <- "data/controls/VDJ_Control_Samples"
files <-  list.files(path = unfiltered.path, pattern = glob2rx("*.adap.txt"))
filtered.files <- list.files(path = filtered.path, pattern = glob2rx("*.tsv"))
controls.files <- list.files(path = controls.path, pattern = glob2rx("*.tsv"))
## Function definitions
##----------------------------------------------------------------------------

## Utility function. Given a filename of the form 'GTSP....[a|b|c]', the function returns the
## GTSP accession number and replicate (where a=1, b=2, c=3)
splitAccnRepl <- function(fname) {
  r <- substr(fname, 9, 9)
  # currently cannot handle replicates past "c". if more, make better soltn
  stopifnot(r %in% c("a", "b", "c"))
  replicate <- ifelse(r == "a", 1,
                      ifelse(r == "b", 2, 3))
  accn <- substr(fname, 1, 8)
  return(c(accn, replicate))
}

## Retrieves information from a MySQL database about a particular sample
## determined by the filename (of the form 'GTSP....[a|b|c]')
##
## Note: functions beginning with '.' are internal functions and
##       are generally run in a rbind/lapply wrapper.
.createMetadata <- function(fname, con) {
  path <- file.path(unfiltered.path, fname)
  accn.repl <- splitAccnRepl(fname)
  accn <- accn.repl[1]
  replicate <- accn.repl[2]
  mdat <- dbGetQuery(con, sprintf("select Patient, CellType, Timepoint from %s where SpecimenAccNum like '%s';",
                                  my.table, accn))
  patient <- mdat$Patient
  cell.type <- mdat$CellType
  timepoint <- mdat$Timepoint
  
  return(data.frame(path=path, accn=accn, patient=patient, cell.type=cell.type,
                    timepoint=timepoint, replicate=replicate))
}

## Reads in unfiltered, whitespace-delimited sequence file with {seq, copy}
## column order and incorporates metadata (from the metadata data frame) into
## the returned data frame.
##
## Note: functions beginning with '.' are internal functions and
##       are generally run in a rbind/lapply wrapper.
.readSeqData <- function(fname, metadata) {
  dat.tmp <- read.delim(as.character(fname), sep="", header=F, 
                        stringsAsFactors=F)
  meta    <- as.list(subset(metadata, path == fname))
  return(data.frame(seq=dat.tmp[,1], copy=dat.tmp[,2], meta))  
}


## Reads in a processed, filtered sequence data from Adaptive
##
## Note: functions beginning with '.' are internal functions and
##       are generally run in a rbind/lapply wrapper.
.readFiltSeqData <- function(fname, path) {
  accn.repl <- splitAccnRepl(fname)
  .accn <- accn.repl[1]
  .replicate <- accn.repl[2]
  dat.tmp <- read.delim(file.path(path, as.character(fname)), stringsAsFactors=F, na.strings="(undefined)")
  # Renaming 'copy' column to avoid collision with unfiltered sequence df column 
  dat.tmp$filt.copy <- dat.tmp$copy
  dat.tmp$seq <- dat.tmp$nucleotide
  dat.tmp$nucleotide <- NULL
  dat.tmp$copy <- NULL
  dat.tmp <- data.frame(accn = .accn, replicate = .replicate, dat.tmp)
  dat.tmp$accn <- as.character(dat.tmp$accn)
  dat.tmp$replicate <- as.character(dat.tmp$replicate)
  return(dat.tmp)
}

.readUnfiltered <- function(metadata) {
  cat("Reading unfiltered sequence data... \n")
  df <- do.call(rbind, lapply(metadata$path, .readSeqData, metadata))
  # Store all the sequences as their reverse complement 
  # (to allow comparison w/ filtered data)
  df$seq <- DNAStringSet(as.character(df$seq))
  df$seq <- as.character(reverseComplement(df$seq))
  # Convert from factors
  df$accn <- as.character(df$accn)
  df$replicate <- as.character(df$replicate)
  return(df)
}

## Write unique sequences to FASTA format from data frame
writeSeqToFasta <- function(df, outname, path=root) {
  .unique <- unique(df$seq)
  names(.unique) <- 1:length(.unique)
  outfile <- file.path(root, outname)
  write.fasta(as.list(.unique), 1:length(.unique), file.out=outfile)
}


readPSLFile <- function(fname) {
  psl <- read.delim(fname, stringsAsFactors=F, header=F)
  colnames(psl) <- c("match", "mismatch", "rep.match", "Ns", "Q.gap.count", 
                     "Q.gap.bases", "T.gap.count", "T.gap.bases", "strand", 
                     "Q.name", "Q.size", "Q.start", "Q.end", "T.name", 
                     "T.size", "T.start", "T.end", "block.count", 
                     "block.sizes", "q.starts", "t.starts")
  len <- psl$Q.end - psl$Q.start
  psl$pct.id <- psl$match / len
  return(psl)
}


idSingletons <- function(df) {
  singletons <- subset(df, copy == 1, select=c(seq, accn))
  across.replicates <- which(duplicated(singletons))
  repl.singletons <- singletons[across.replicates, ]
  n.singletons.replicated <- length(across.replicates)
  n.singletons <- length(singletons$seq)
  proportion.singletons <- n.singletons.replicated/n.singletons
  cat("\n\t", "Total singletons across replicates:", n.singletons.replicated,
      "\n\t", "Proportion of total singletons:", n.singletons.replicated, "/", 
      n.singletons, "=", proportion.singletons, "\n\n")
  print(summary(repl.singletons))
}

## Generates the master data frame with the aggregated sequence info. 
if (TRUE) {

  cat("Collecting metadata... \n")
  con <- dbConnect(MySQL(), host=my.host, user=my.user, password=my.pass, 
                   dbname=my.db)
  metadata  <- do.call(rbind, lapply(files, .createMetadata, con))

  cat("Reading unfiltered sequence data... \n")
  df.unfilt <- do.call(rbind, lapply(metadata$path, .readSeqData, metadata))
  # Store all the sequences as their reverse complement 
  # (to allow comparison w/ filtered data)
  df.unfilt$seq <- DNAStringSet(as.character(df.unfilt$seq))
  df.unfilt$seq <- as.character(reverseComplement(df.unfilt$seq))
  # Convert from factors
  df.unfilt$accn <- as.character(df.unfilt$accn)
  df.unfilt$replicate <- as.character(df.unfilt$replicate)

  cat("Reading filtered sequence data... \n")
  df.filt <- do.call(rbind, lapply(filtered.files, .readFiltSeqData, 
                                   filtered.path))
  
  cat("Reading control sequence data... \n")
  df.ctrl <- do.call(rbind, lapply(controls.files, .readFiltSeqData,
                                   controls.path))

  cat("Joining data by sequence, replicate and accession... \n")
  df.filt.all <- merge(df.filt, df.ctrl, all=T)
  ap <- df.filt.all
  save(ap, file="all_processed.RData")
  # change this to TRUE to continue with the rest of the processing (obsolete)
  stopifnot(F)
  df <- merge(df.unfilt, df.filt.all, by=c("accn", "seq", "replicate"), all=T)
  # idx is an unique index for each row that is used to identify that row
  # when the sequences are processed by external scripts 
  df$idx <- 1:length(df$seq)
  cat("Writing unique sequences as .fasta file...\n")
  unique.seq <- unique(df$seq)
  fasta.file = 'unique.fa'
  writeSeqToFasta(df, fasta.file)
  # adds index to each sequence to re-merge with BLAT results
  unique.seq <- data.frame(cbind(unique.seq, 1:length(unique.seq)))
  colnames(unique.seq) = c("seq", "idx")
  unique.seq$idx <- as.integer(as.numeric(unique.seq$idx))
  # This takes upwards of a couple hours. If it's already done, just ensure
  # that the variable `psl.file` points to the correct psl output.
  if (FALSE) {
    cat("BLAT'ing all sequences against each other...\n")
    cat("This may take a very long time. Grab a Snickers.\n")
    min.score = 50
    system(sprintf('python blat_runner.py -m %i %s', min.score, fasta.file))
  } 

  psl.file <- 'unique.psl'

  cat("Processing psl file from", psl.file, '...\n')
  processed.file <- 'unique.processed'
  psl.max = 0.99
  psl.min = 0.84
  system(sprintf('python psl_processing.py --max %1.2f --min %1.2f %s', 
                 psl.max, psl.min, psl.file))

  cat("Reading in processed psl file (see psl_processing.py)...")
  psl <- read.delim(processed.file, stringsAsFactors=F)
  psl$idx <- as.integer(as.numeric(psl$idx))
  psl <- merge(psl, unique.seq, by=c("idx"))
  psl$idx <- NULL
  df <- merge(df, psl, by=c("seq"))
  idcols <- colnames(psl)[!colnames(psl) %in% c("seq", "idx")]
  
  cat("done.\n")
  
  all.data <- df
  filename <- 'all_data.RData'
  cat("Saved file as", filename)
  save(all.data, file = filename)

  cat("Starting analysis:")
  
  .pctAboveThreshold <- function(threshold, .df, gt.matches=1) {
    l = length(.df$idx)
    a = length(.df$idx[.df[, threshold] > gt.matches])
    return(a/l)
  }
  
  # Returns the percentage of sequences in the given copy number range
  # that match other sequences at the thresholds given in the psl_processing.py
  # script (by default 99-84% identity at 5% intervals)
  matchesAboveCopyNumbers <- function(df, cn.range) {
    matches <- data.frame(sapply(cn.range, 
                                 function(x) sapply(idcols, 
                                                    .pctAboveThreshold, 
                                                    subset(df, copy==x))))
    colnames(matches) <- c(paste("cn", cn.range, sep='.'))
    return(matches)
  }
  
  matches <- matchesAboveCopyNumbers(df, 1:4)

}






