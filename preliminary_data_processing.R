##----------------------------------------------------------------------------
##  Data processing functions for unfiltered VDJ sequencing data
##  Erik Clarke (ecl@mail.med.upenn.edu)
##----------------------------------------------------------------------------
source("misc_functions.R")
path  <- 	"data/patients/processed/unfiltered"
files <-  list.files(path = path, pattern = glob2rx("*.adap.txt"))


## Function definitions
##----------------------------------------------------------------------------
createMetadata <- function(string) {
  r <- substr(string, 9, 9)
  stopifnot(r %in% c("a", "b", "c"))
  replicate <- ifelse(r == "a", 1,
                      ifelse(r == "b", 2, 3))
  patient <- substr(string, 1, 8)
  path <- file.path(path, string)
  accn <- paste(substr(string, 1, 8), replicate, sep='.')
  return(data.frame(path=path, patient=patient, replicate=replicate, accn=accn))
}

readSeqData <- function(file) {
  dat.tmp <- read.delim(as.character(file), sep="", header=F, 
                        stringsAsFactors=F)
  accn    <- metadata$accn[metadata$path == file]
  return(data.frame(seq=dat.tmp[,1], copy=dat.tmp[,2], accn=accn))  
}


## Script body
##----------------------------------------------------------------------------
## Data importing
cat("Collecting paths to data...\n")
metadata  <- do.call(rbind, lapply(files, createMetadata))
cat("Reading sequence data... ")
data      <- do.call(rbind, lapply(metadata$path, readSeqData))
data$seq  <- as.character(df$seq)
cat("done.")

## Identifying singletons


