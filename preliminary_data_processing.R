##----------------------------------------------------------------------------
##  Data processing functions for unfiltered VDJ sequencing data
##  Erik Clarke (ecl@mail.med.upenn.edu)
##----------------------------------------------------------------------------
library(RMySQL)
library(seqinr)


## Constants
##----------------------------------------------------------------------------

# Load MySQL connection info
# This should have variables named my.user, my.pass, my.host, my.db, my.table
source("mysql_info.R")

# File info
path  <- 	"data/patients/processed/unfiltered"
files <-  list.files(path = path, pattern = glob2rx("*.adap.txt"))


## Function definitions
##----------------------------------------------------------------------------

createMetadata <- function(string, con) {
  path <- file.path(path, string)
  r <- substr(string, 9, 9)
  # currently cannot handle replicates past "c". if more, make better soltn
  stopifnot(r %in% c("a", "b", "c"))
  replicate <- ifelse(r == "a", 1,
                      ifelse(r == "b", 2, 3))
  accn <- substr(string, 1, 8)
  mdat <- dbGetQuery(con, sprintf("select Patient, CellType, Timepoint from %s where SpecimenAccNum like '%s';",
                                  my.table, accn))
  patient <- mdat$Patient
  cell.type <- mdat$CellType
  timepoint <- mdat$Timepoint
  
  return(data.frame(path=path, accn=accn, patient=patient, cell.type=cell.type,
                    timepoint=timepoint, replicate=replicate))
}

# Reads in unfiltered, whitespace-delimited sequence file with {seq, copy}
# column order (use with lapply, as below)
readSeqData <- function(file) {
  dat.tmp <- read.delim(as.character(file), sep="", header=F, 
                        stringsAsFactors=F)
  meta    <- as.list(subset(metadata, path == file))
  return(data.frame(seq=dat.tmp[,1], copy=dat.tmp[,2], meta))  
}


## Script
##----------------------------------------------------------------------------

## Read sequence data
cat("Collecting metadata... ")
con <- dbConnect(MySQL(), host=my.host, user=my.user, password=my.pass, dbname=my.db)
metadata  <- do.call(rbind, lapply(files, createMetadata, con))
cat("done.\n")
cat("Reading sequence data... ")
df <- do.call(rbind, lapply(metadata$path, readSeqData))
df$seq <- as.character(df$seq)
cat("done.\n")


## Write sequence data in FASTA format
## This takes a while and only needs to be done once.
if (FALSE) {
  seqnames <- paste(df$accn, 1:length(df$seq), sep='_')
  seqs <- df$seq
  names(seqs) <- seqnames
  fasta.file <- file.path(path, "unfiltered.all.fasta")
  write.fasta(as.list(seqs), seqnames, file.out=fasta.file)
}

## Run igblast all the sequences against the immunoglobulin V/D/J genes.
## (This takes considerable time, so ensure you want to do it by
##  switching condition to TRUE)
if (FALSE) {
  # Ensure you have igblast installed from NCBI with the database directory
  # included (or specify path below)
  dbpath <- "/usr/local/ncbi/igblast/database/"
  fasta.file <- file.path(path, "unfiltered.all.fasta")
  out.file <- file.path(path, "unfiltered.all.blast.out")
  cmd <- sprintf("igblastn -germline_db_V %s -germline_db_D %s -germline_db_J %s -organism human -domain_system kabat -query %s -out %s",
                 file.path(dbpath, "human_gl_V"), file.path(dbpath, "human_gl_D"), file.path(dbpath, "human_gl_J", fasta.file), out.file)
  system(cmd)
}



## Identify singletons appearing across replicates
singletons <- subset(df, copy == 1, select=c(seq, accn))
across.replicates <- which(duplicated(singletons))
repl.singletons <- singletons[across.replicates, ]
# Output results
n.singletons.replicated <- length(across.replicates)
n.singletons <- length(singletons$seq)
proportion.singletons <- n.singletons.replicated/n.singletons
cat("\n\t", "Total singletons across replicates:", n.singletons.replicated,
    "\n\t", "Proportion of total singletons:", n.singletons.replicated, "/", 
    n.singletons, "=", proportion.singletons, "\n\n")
print(summary(repl.singletons))

## Identify single- or doubletons (lowcn) appearing across cell types
# Limit to samples at the d365 timepoint
tmp.df <- subset(df, timepoint == "d365")
for (.accn in unique(tmp.df$accn)){
  df.rep1 <- subset(tmp.df, replicate == 1 && accn == .accn)
  for (r in 2:3) {
    sset <- subset(tmp.df, replicate == r && accn == .accn)
    for (.seq in sset$seq) {
      row <- subset(sset, seq == .seq)
      df.rep1[df.rep1$seq == .seq]$copy += row$copy
    }
  }
}
df.rep1 <- subset(tmp.df, replicate == 1)
for (r in 2:3) {
  sset <- subset(tmp.df, replicate == r)
  for (seq in sset$seq) {
    row <- subset(sset, seq == seq)
    df.rep1[df.rep1$seq == seq]$copy
  }
}
for (seq in tmp.df$seq) {
  r <- subset(tmp.df, seq == seq)
  if (r$replicate > 1) {
    copy <-
  }
} 
singletons <- subset(df, copy <= 2, select = c(seq, accn, cell.type))
#lowcn <- subset(df, copy <= 2, select=c())







