##----------------------------------------------------------------------------
##  Data processing functions for unfiltered VDJ sequencing data
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
root <- "/Volumes/ecl/R/VJD_XSCID"
unfiltered.path  <-   "data/patients/processed/unfiltered"
filtered.path <- "data/patients/processed/filtered"
files <-  list.files(path = unfiltered.path, pattern = glob2rx("*.adap.txt"))
filtered.files <- list.files(path = filtered.path, pattern = glob2rx("*.tsv"))

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
  dat.tmp <- read.delim(file.path(path, as.character(fname)), stringsAsFactors=F)
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

## Write sequence data in FASTA format
writeSeqToFasta <- function(df, outname, path=root) {
  seqnames <- paste(df$accn, df$replicate, df$idx, sep='_')
  seqs <- df$seq
  names(seqs) <- seqnames
  fasta.file <- file.path(root, outname)
  write.fasta(as.list(seqs), seqnames, file.out=fasta.file)
}

## Identify singletons that
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
  con <- dbConnect(MySQL(), host=my.host, user=my.user, password=my.pass, dbname=my.db)
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
  df.filt <- do.call(rbind, lapply(filtered.files, .readFiltSeqData, filtered.path))

  cat("Joining data by sequence, replicate and accession... \n")
  df <- merge(df.unfilt, df.filt, by=c("accn", "seq", "replicate"), all=T)
  # idx is an unique identifier for each row that is used to identify that row
  # when the sequences are processed by external programs (i.e. as a FASTA accn #)
  df$idx <- 1:length(df$seq)
  
  cat("Reading in processed psl file (see psl_processing.py)...")
  psl <- read.delim(file.path(root, "ig.all.psl.processed"))
  df <- merge(df, psl, by=c("accn", "idx"))
  idcols <- colnames(psl)[!colnames(psl) %in% c("accn", "idx")]
  
  .pctAboveThreshold <- function(threshold, .df, gt.matches=1) {
    l = length(.df$idx)
    a = length(.df$idx[.df[, threshold] > gt.matches])
    cat("a/l:", a, "/", l, "\t")
    return(a/l)
  }
  
  matchesAboveCopyNumbers <- function(df, cn.range = 1:3) {
    matches <- data.frame(sapply(cn.range, 
                                 function(x) sapply(idcols, .pctAboveThreshold, subset(df, copy==x))))
    colnames(matches) <- c(paste("cn", cn.range, sep='.'))
    return(matches)
  }
  
  print(matchesAboveCopyNumbers(df, 1:3))
  
  matches <- data.frame(sapply(1:3, function(x) 
    sapply(idcols, .pctAboveThreshold, subset(df, copy==x))))
  colnames(matches) <- c(paste("cn", 1:3, sep='.'))
  for (col in idcols) {
    cat(col)
    threshold = substr(col, 4, 5)
    for (cn in 1:3){
      .tmp = subset(df, copy==cn)
      cat(dim(.tmp))
      cat("\nSequences with copy number", cn, ":", length(.tmp$idx))
      cat("\n\tNumber of these with > 1 matches at", threshold, "%:", length(.tmp$idx[.tmp[,col] > 1]))
      cat("\n")
    }
  }
  cat("done.\n")
}






